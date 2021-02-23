suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing id and species columns"),

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
counts_file <- opt$counts_file
alignments_file <- opt$alignments_file
adapter_alignments_file <- opt$adapter_alignments_file
output_prefix <- opt$output_prefix

if (is.null(samples_file)) stop("Samples file must be specified")
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
missing_columns <- setdiff(c("id", "fastq", "species", "genome", "control genome", "id_prefix"), colnames(samples))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", samples_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
samples <- mutate(samples, id = as_factor(id))
if (nrow(filter(samples, is.na(samples$id))) > 0) {
  stop("Missing sample identifiers")
}
if (nlevels(samples$id) != nrow(samples)) {
  stop("Duplicate sample identifiers")
}

# read sequence counts summary from FASTQ sampling
message("Reading sequence counts from sampling")
counts <- read_csv(counts_file, col_types = "cii")
missing_columns <- setdiff(c("id", "read", "sampled"), colnames(counts))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", counts_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
counts <- select(counts, id_prefix = id, sequences = read, sampled)
samples <- left_join(samples, counts, by = "id_prefix")
missing_counts <- samples %>%
  select(id, sequences, sampled) %>%
  filter(!complete.cases(.))
if (nrow(missing_counts) > 0) {
  stop("Missing counts for sample(s): ", str_c(missing_counts$id, collapse = ", "))
}

# output sequence counts file
samples %>%
  select(id, fastq, species, genome, control, `control genome`, sequences, sampled) %>%
  write_csv(output_counts_file, na = "")

# read bowtie alignment output
message("Reading genome alignments")
alignments <- read_tsv(alignments_file, col_types = "cccciccic")
alignment_columns <- c("genome", "read", "strand", "chromosome", "position", "sequence", "quality", "num", "mismatches")
missing_columns <- setdiff(alignment_columns, colnames(alignments))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# calculate sequence lengths and count mismatches
message("Computing sequence lengths and mismatch/error rates")
alignments <- alignments %>%
  mutate(length = nchar(sequence)) %>%
  mutate(mismatch_count = str_count(mismatches, ",") + 1L) %>%
  mutate(mismatch_count = as.integer(replace_na(mismatch_count, 0L)))

# separate sample name from read id
message("Separating read identifiers and matching to samples")
alignments <- alignments %>%
  separate(read, into = c("id_prefix", "read"), sep = "\\|", extra = "merge")

alignment_count <- nrow(alignments)

alignments <- alignments %>%
  inner_join(select(samples, id, id_prefix), by = "id_prefix") %>%
  select(id, read, genome, chromosome, position, strand, sequence, quality, length, mismatches = mismatch_count)

if (nrow(alignments) != alignment_count) {
  stop("Found alignments that couldn't be matched to a sample")
}

# mark alignments to the expected genome(s)
message("Determining matches to expected and control genomes")
expected_genomes <- samples %>%
  select(id, genome) %>%
  filter(!is.na(genome)) %>%
  mutate(expected = TRUE)
control_genomes <- samples %>%
  select(id, genome = `control genome`) %>%
  filter(!is.na(genome)) %>%
  mutate(control = TRUE)
alignments <- alignments %>%
  left_join(expected_genomes, by = c("id", "genome")) %>%
  left_join(control_genomes, by = c("id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  mutate(control = replace_na(control, FALSE))

# alignment summary
aligned_summary <- alignments %>%
  group_by(id, genome, expected, control) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(samples, id, total = sampled), by = "id") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(id, genome, expected, control, count, percentage, error_rate)

# mark alignments with fewest mismatches
message("Determining best genomic alignments for each read")
if (nrow(alignments) == 0) {
  alignments <- mutate(alignments, fewest_mismatches = FALSE)
} else {
  alignments <- alignments %>%
    group_by(id, read) %>%
    mutate(fewest_mismatches = mismatches == min(mismatches)) %>%
    ungroup()
}

fewest_mismatch_counts <- alignments %>%
  filter(fewest_mismatches) %>%
  count(id, genome, name = "fewest_mismatch_count")

# assign reads to genomes/species
message("Assigning reads to genomes")
alignments <- alignments %>%
  left_join(fewest_mismatch_counts, by = c("id", "genome")) %>%
  group_by(id, read) %>%
  arrange(desc(expected), desc(fewest_mismatch_count)) %>%
  mutate(assigned = row_number() == 1) %>%
  ungroup() %>%
  select(-fewest_mismatches, -fewest_mismatch_count)

assigned_alignments <- filter(alignments, assigned)

assigned_summary <- assigned_alignments %>%
  group_by(id, genome, expected, control) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(samples, id, total = sampled), by = "id") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(id, genome, expected, control, count, percentage, error_rate)

# output alignments
alignments %>%
  select(-length) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_tsv(output_alignments_file, na = "")

# unmapped read counts
unmapped_counts <- assigned_alignments %>%
  count(id, name = "aligned") %>%
  full_join(select(samples, id, sampled), by = "id") %>%
  mutate(aligned = replace_na(aligned, 0)) %>%
  mutate(count = sampled - aligned) %>%
  mutate(percentage = ifelse(sampled == 0, 0.0, round(100 * count / sampled, digits = 1))) %>%
  select(id, count, percentage)

# adapter alignments
message("Reading adapter alignments")
adapter_alignments <- read_tsv(adapter_alignments_file, col_types = "ciicciicdi")
adapter_alignment_columns <- c("read", "start", "end", "strand", "adapter", "adapter start", "adapter end", "adapter strand", "percent identity", "score")
missing_columns <- setdiff(adapter_alignment_columns, colnames(adapter_alignments))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", adapter_alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

adapter_alignment_count <- nrow(adapter_alignments)

adapter_alignments <- adapter_alignments %>%
  separate(read, into = c("id_prefix", "read"), sep = "\\|", extra = "merge") %>%
  inner_join(select(samples, id, id_prefix), by = "id_prefix") %>%
  select(id, everything(), -id_prefix)

if (nrow(adapter_alignments) != adapter_alignment_count) {
  stop("Found adapter alignments that couldn't be matched to a sample")
}

write_tsv(adapter_alignments, output_adapter_alignments_file, na = "")

adapter_counts <- adapter_alignments %>%
  count(id, name = "count") %>%
  full_join(select(samples, id, sampled), by = "id") %>%
  mutate(count = replace_na(count, 0)) %>%
  mutate(percentage = ifelse(sampled == 0, 0.0, round(100 * count / sampled, digits = 1))) %>%
  select(id, count, percentage)

# summary
message("Creating alignment summary table")
summary <- aligned_summary %>%
  select(id, genome, expected, control, aligned = count, `aligned %` = percentage, `aligned error rate` = error_rate) %>%
  left_join(select(assigned_summary, id, genome, assigned = count, `assigned %` = percentage, `assigned error rate` = error_rate), by = c("id", "genome")) %>%
  bind_rows(transmute(unmapped_counts, id, genome = "unmapped", aligned = count, `aligned %` = percentage, assigned = count, `assigned %` = percentage)) %>%
  mutate(assigned = replace_na(assigned, 0L)) %>%
  mutate(`assigned %` = replace_na(`assigned %`, 0.0)) %>%
  bind_rows(transmute(adapter_counts, id, genome = "adapters", aligned = count, `aligned %` = percentage)) %>%
  mutate(id = as_factor(id))

summary %>%
  arrange(id, is.na(expected), desc(assigned)) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_csv(output_summary_file, na = "")

# stacked bar plots
message("Creating bar charts")

pdf(output_plot_file, width = 8, height = 1 + nrow(distinct(summary, id)) * 0.55)

lower_error_rate <- 0.25
upper_error_rate <- 1.5

summary %>%
  left_join(select(samples, id, sequences, sampled), by = "id") %>%
  filter(!is.na(assigned)) %>%
  filter(assigned > 0) %>%
  mutate(category = ifelse(expected, "expected species", "unexpected/contaminant")) %>%
  mutate(category = ifelse(control, "control", category)) %>%
  mutate(category = replace_na(category, "unmapped")) %>%
  mutate(category = factor(category, levels = rev(c("expected species", "control", "unexpected/contaminant", "unmapped")))) %>%
  mutate(id = fct_rev(id)) %>%
  mutate(sequences = (assigned / sampled) * sequences / 1e6) %>%
  mutate(bounded_error_rate = pmin(`assigned error rate`, upper_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = pmax(bounded_error_rate, lower_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = ifelse(category == "unmapped", lower_error_rate, bounded_error_rate)) %>%
  arrange(id, category, desc(sequences)) %>%
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
  mutate(adapter = count, `no adapter` = sampled - count) %>%
  mutate_at(vars(adapter, `no adapter`), ~ (. / sampled) * sequences / 1e6) %>%
  select(id, adapter, `no adapter`) %>%
  pivot_longer(-id, names_to = "category", values_to = "sequences") %>%
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


