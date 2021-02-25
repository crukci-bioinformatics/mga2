suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing details of species and controls for each sample/dataset"),
  
  make_option(c("--genomes"), dest = "genomes_file",
              help = "Genomes file containing list of indexed reference genome sequences and associated species"),
  
  make_option(c("--counts"), dest = "counts_file",
              help = "Sequence counts file output by sample-fastq containing id, read and sampled columns"),

  make_option(c("--alignments"), dest = "alignments_file",
              help = "Alignments file with collated output from bowtie with genome, read, strand, chromosome, position, sequence, quality, num and mismatch columns"),
  
  make_option(c("--adapter-alignments"), dest = "adapter_alignments_file",
              help = "Adapter alignments file with collated output from exonerate with read, start, end, strand, adapter, adapter_start, adapter_end, adapter_strand, percent_identity and score columns"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

samples_file <- opt$samples_file
genomes_file <- opt$genomes_file
counts_file <- opt$counts_file
alignments_file <- opt$alignments_file
adapter_alignments_file <- opt$adapter_alignments_file
output_prefix <- opt$output_prefix

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(genomes_file)) stop("Genomes file must be specified")
if (is.null(counts_file)) stop("Sequence counts file must be specified")
if (is.null(alignments_file)) stop("Alignments file must be specified")
if (is.null(adapter_alignments_file)) stop("Adapter alignments file must be specified")

if (is.null(output_prefix)) output_prefix <- "mga."

suppressPackageStartupMessages(library(tidyverse))

output_counts_file <- str_c(output_prefix, "sequence_counts.csv")
output_alignments_file <- str_c(output_prefix, "alignments.txt")
output_adapter_alignments_file <- str_c(output_prefix, "adapter_alignments.txt")
output_summary_file <- str_c(output_prefix, "summary.csv")
output_plot_file <- str_c(output_prefix, "bar_charts.pdf")

# read sample sheet
message("Reading sample information")
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))

expected_columns <- c("rownum", "id", "fastq", "species", "controls", "genomes", "control_genomes")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns in ", samples_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
samples <- select(samples, all_of(expected_columns))

if (nrow(filter(samples, is.na(rownum))) > 0 || nrow(filter(samples, is.na(id))) > 0) {
  stop("missing sample identifiers")
}

duplicates <- samples %>%
  count(rownum) %>%
  filter(n > 1) %>%
  pull(rownum)
if (length(duplicates) > 0) {
  stop("duplicate id prefixes found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

duplicates <- samples %>%
  count(id) %>%
  filter(n > 1) %>%
  pull(id)
if (length(duplicates) > 0) {
  stop("duplicate ids found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# read genomes that have been aligned to
message("Reading genome details")
genomes <- read_csv(genomes_file, col_types = cols(.default = col_character()))

expected_columns <- c("genome", "species")
missing_columns <- setdiff(expected_columns, colnames(genomes))
if (length(missing_columns) > 0) {
  stop("missing columns in ", genomes_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
genomes <- select(genomes, all_of(expected_columns))

# read sequence counts summary from FASTQ sampling
message("Reading sequence counts from sampling")
counts <- read_csv(counts_file, col_types = "cii")

expected_columns <- c("id", "read", "sampled")
missing_columns <- setdiff(expected_columns, colnames(counts))
if (length(missing_columns) > 0) {
  stop("missing columns in ", counts_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
counts <- select(counts, rownum = id, sequences = read, sampled)

samples <- left_join(samples, counts, by = "rownum")
missing_counts <- samples %>%
  select(id, sequences, sampled) %>%
  filter(!complete.cases(.))
if (nrow(missing_counts) > 0) {
  stop("missing counts for sample(s): ", str_c(missing_counts$id, collapse = ", "))
}

# output sequence counts file
samples %>%
  select(id, fastq, species, genomes, controls, `control genomes` = control_genomes, sequences, sampled) %>%
  write_csv(output_counts_file, na = "")

# read bowtie alignment output
message("Reading genome alignments")
alignments <- read_tsv(alignments_file, col_types = "cccciccic")
expected_columns <- c("genome", "read", "strand", "chromosome", "position", "sequence", "quality", "num", "mismatches")
missing_columns <- setdiff(expected_columns, colnames(alignments))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
alignments <- select(alignments, all_of(expected_columns))

# separate sample id prefix from read id
message("Separating read identifiers and matching to samples")
alignment_count <- nrow(alignments)
alignments <- alignments %>%
  separate(read, into = c("rownum", "read"), sep = "\\|", extra = "merge") %>%
  inner_join(select(samples, id, rownum), by = "rownum") %>%
  select(rownum, id, read, genome, chromosome, position, strand, sequence, quality, mismatches)
if (nrow(alignments) != alignment_count) {
  stop("found alignments that couldn't be matched to a sample")
}

# calculate sequence lengths and count mismatches
message("Computing sequence lengths and mismatch/error rates")
alignments <- alignments %>%
  mutate(length = nchar(sequence)) %>%
  mutate(mismatches = str_count(mismatches, ",") + 1L) %>%
  mutate(mismatches = as.integer(replace_na(mismatches, 0L)))

# mark alignments to the expected genome(s)
message("Determining matches to expected and control genomes")

expected_genomes <- samples %>%
  select(rownum, id, genome = genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na()

control_genomes <- samples %>%
  select(rownum, id, genome = control_genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na()

expected_genomes <- expected_genomes %>%
  bind_rows(control_genomes) %>%
  mutate(expected = TRUE)

control_genomes <- mutate(control_genomes, control = TRUE)

alignments <- alignments %>%
  left_join(expected_genomes, by = c("rownum", "id", "genome")) %>%
  left_join(control_genomes, by = c("rownum", "id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  mutate(control = replace_na(control, FALSE))

# mark alignments with fewest mismatches
message("Determining best genomic alignments for each read")
if (nrow(alignments) == 0) {
  alignments <- mutate(alignments, fewest_mismatches = FALSE)
} else {
  alignments <- alignments %>%
    group_by(rownum, read) %>%
    mutate(fewest_mismatches = mismatches == min(mismatches)) %>%
    ungroup()
}

fewest_mismatch_counts <- alignments %>%
  filter(fewest_mismatches) %>%
  count(rownum, genome, name = "fewest_mismatch_count")

# assign reads to genomes/species
message("Assigning reads to genomes")
alignments <- alignments %>%
  left_join(fewest_mismatch_counts, by = c("rownum", "genome")) %>%
  group_by(rownum, id, read) %>%
  arrange(desc(expected), desc(fewest_mismatch_count)) %>%
  mutate(assigned = row_number() == 1) %>%
  ungroup() %>%
  select(-fewest_mismatches, -fewest_mismatch_count)

# output alignments
alignments %>%
  arrange(rownum, id, read, genome) %>%
  select(-rownum, -length) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_tsv(output_alignments_file, na = "")

# adapter alignments
message("Reading adapter alignments")
adapter_alignments <- read_tsv(adapter_alignments_file, col_types = "ciicciicdi")
expected_columns <- c("read", "start", "end", "strand", "adapter", "adapter start", "adapter end", "adapter strand", "percent identity", "score")
missing_columns <- setdiff(expected_columns, colnames(adapter_alignments))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", adapter_alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
adapter_alignments <- select(adapter_alignments, all_of(expected_columns))

# separate sample id prefix from read id
message("Separating read identifiers for adapter alignments and matching to samples")
adapter_alignment_count <- nrow(adapter_alignments)
adapter_alignments <- adapter_alignments %>%
  separate(read, into = c("rownum", "read"), sep = "\\|", extra = "merge") %>%
  inner_join(select(samples, id, rownum), by = "rownum") %>%
  select(rownum, id, everything())
if (nrow(adapter_alignments) != adapter_alignment_count) {
  stop("found adapter alignments that couldn't be matched to a sample")
}

# output adapter alignments
adapter_alignments %>%
  arrange(rownum, id, read, start) %>%
  select(-rownum) %>%
  write_tsv(output_adapter_alignments_file, na = "")

# alignment summary
message("Creating alignment summary table")

aligned_summary <- alignments %>%
  group_by(rownum, genome) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(samples, rownum, id, total = sampled), by = "rownum") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(rownum, id, genome, count, percentage, error_rate) %>%
  full_join(expected_genomes, by = c("rownum", "id", "genome")) %>%
  full_join(control_genomes, by = c("rownum", "id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  mutate(control = replace_na(control, FALSE)) %>%
  select(rownum, id, genome, expected, control, count, percentage, error_rate)

assigned_alignments <- filter(alignments, assigned)

assigned_summary <- assigned_alignments %>%
  group_by(rownum, genome) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(samples, rownum, total = sampled), by = "rownum") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(rownum, genome, count, percentage, error_rate)

# unmapped read counts
unmapped_counts <- assigned_alignments %>%
  count(rownum, name = "aligned") %>%
  full_join(select(samples, rownum, id, sampled), by = "rownum") %>%
  mutate(aligned = replace_na(aligned, 0)) %>%
  mutate(count = sampled - aligned) %>%
  mutate(percentage = ifelse(sampled == 0, 0.0, round(100 * count / sampled, digits = 1))) %>%
  select(rownum, id, count, percentage)

# counts of reads that align to adapters
adapter_counts <- adapter_alignments %>%
  count(rownum, name = "count") %>%
  full_join(select(samples, rownum, id, sampled), by = "rownum") %>%
  mutate(count = replace_na(count, 0)) %>%
  mutate(percentage = ifelse(sampled == 0, 0.0, round(100 * count / sampled, digits = 1))) %>%
  select(rownum, id, count, percentage)

# combined summary table

summary <- aligned_summary %>%
  left_join(genomes, by = "genome") %>%
  select(rownum, id, genome, species, expected, control, aligned = count, `aligned %` = percentage, `aligned error rate` = error_rate) %>%
  left_join(select(assigned_summary, rownum, genome, assigned = count, `assigned %` = percentage, `assigned error rate` = error_rate), by = c("rownum", "genome")) %>%
  bind_rows(transmute(unmapped_counts, rownum, id, genome = "unmapped", aligned = count, `aligned %` = percentage, assigned = count, `assigned %` = percentage)) %>%
  mutate(across(c("aligned", "aligned %", "assigned", "assigned %"), ~ replace_na(., 0L))) %>%
  bind_rows(transmute(adapter_counts, rownum, id, genome = "adapters", aligned = count, `aligned %` = percentage))

summary %>%
  arrange(rownum, is.na(expected), desc(assigned)) %>%
  select(-rownum) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_csv(output_summary_file, na = "")

# stacked bar plots
message("Creating bar charts")

pdf(output_plot_file, width = 8, height = 1 + nrow(distinct(summary, id)) * 0.55)

lower_error_rate <- 0.25
upper_error_rate <- 1.5

summary %>%
  filter(!is.na(assigned)) %>%
  left_join(select(samples, rownum, sequences, sampled), by = "rownum") %>%
  mutate(sequences = ifelse(sequences == 0, 0, (assigned / sampled) * sequences / 1e6)) %>%
  mutate(category = ifelse(expected, "expected species", "unexpected/contaminant")) %>%
  mutate(category = ifelse(control, "control", category)) %>%
  mutate(category = replace_na(category, "unmapped")) %>%
  mutate(category = factor(category, levels = rev(c("expected species", "control", "unexpected/contaminant", "unmapped")))) %>%
  arrange(rownum) %>%
  mutate(id = fct_rev(id)) %>%
  mutate(bounded_error_rate = pmin(`assigned error rate`, upper_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = pmax(bounded_error_rate, lower_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = ifelse(category == "unmapped", lower_error_rate, bounded_error_rate)) %>%
  arrange(id, category, sequences) %>%
  mutate(group = row_number()) %>%
  ggplot(aes(x = id, y = sequences, group = group, fill = category, alpha = -bounded_error_rate)) +
  geom_col(colour = "grey30", width = 0.6) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  scale_fill_manual(values = c("grey97", "red", "goldenrod1", "green"), drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_alpha(guide = "none") +
  labs(y = "sequence reads (millions)") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank()
  )

adapter_counts %>%
  left_join(select(samples, id, sequences, sampled), by = "id") %>%
  select(rownum, id, sequences, sampled, adapter = count) %>%
  mutate(`no adapter` = sampled - adapter) %>%
  pivot_longer(all_of(c("adapter", "no adapter")), names_to = "category", values_to = "count") %>%
  mutate(sequences = ifelse(sequences == 0, 0, (count / sampled) * sequences / 1e6)) %>%
  arrange(rownum) %>%
  mutate(id = fct_rev(id)) %>%
  mutate(category = fct_rev(category)) %>%
  ggplot(aes(x = id, y = sequences, fill = category)) +
  geom_col(colour = "grey30", width = 0.6) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  scale_fill_manual(values = c("grey97", "magenta"), drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_alpha(guide = "none") +
  labs(y = "sequence reads (millions)") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank()
  )

dev.off()


