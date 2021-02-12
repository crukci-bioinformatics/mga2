suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing id and species columns"),

  make_option(c("--genomes"), dest = "genomes_file",
              help = "Genomes file containing id and species columns"),
  
  make_option(c("--counts"), dest = "counts_file",
              help = "Sequence counts file output by sample-fastq containing id, read and sampled columns"),

  make_option(c("--alignments"), dest = "alignments_file",
              help = "Alignments file with collated output from bowtie with genome, read, strand, chr, pos, sequence, quality, num and mismatch columns"),

  make_option(c("--control"), dest = "control_species",
              help = "The control or spike-in genome/species (optional)"),

  make_option(c("--summary"), dest = "summary_file",
              help = "Summary CSV output file"),
  
  make_option(c("--plot"), dest = "plot_file",
              help = "Plot PDF output file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

# args <- c(
#   "--samples=samplesheet.csv",
#   "--genomes=genomes.csv",
#   "--counts=sequence_counts.csv",
#   "--alignments=bowtie_alignments.txt",
#   "--control=phix",
#   "--summary=summary.csv",
#   "--plot=summary.pdf"
# )

opt <- parse_args(option_parser, args)

samples_file <- opt$samples_file
genomes_file <- opt$genomes_file
counts_file <- opt$counts_file
alignments_file <- opt$alignments_file
control_species <- opt$control_species
summary_file <- opt$summary_file
plot_file <- opt$plot_file

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(genomes_file)) stop("Genomes file must be specified")
if (is.null(counts_file)) stop("Sequence counts file must be specified")
if (is.null(alignments_file)) stop("Alignments file must be specified")
if (is.null(summary_file)) stop("Summary CSV output file must be specified")
if (is.null(plot_file)) stop("Plot PDF output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
missing_columns <- setdiff(c("id", "species"), colnames(samples))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", samples_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# read genomes file
genomes <- read_csv(genomes_file, col_types = cols(.default = col_character()))
missing_columns <- setdiff(c("genome", "species"), colnames(genomes))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", genomes_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# create synonym lookup table
synonyms <- bind_rows(
  genomes %>%
    select(genome, synonym = genome),
  genomes %>%
    select(genome, synonym = species),
  genomes %>%
    select(genome, synonym = synonyms) %>%
    separate_rows(synonym, sep = "\\|")
) %>%
  drop_na() %>%
  mutate(synonym = str_to_lower(synonym)) %>%
  distinct()

duplicate_synonyms <- synonyms %>%
  count(synonym, name = "count") %>%
  filter(count > 1)

if (nrow(duplicate_synonyms) > 0) {
  duplicate_synonyms %>%
    select(synonym) %>%
    left_join(synonyms, by = "synonym") %>%
    left_join(genomes, by = "genome") %>%
    arrange(synonym, genome) %>%
    print(n = Inf)
  stop("Error: ", genomes_file, " contains the same synonym(s) for more than one genome")
}

# expected species/genomes
expected_species <- samples %>%
  select(id, species) %>%
  mutate(control = FALSE)
if (!is.na(control_species)) {
  expected_species <- bind_rows(
    expected_species,
    mutate(expected_species, species = control_species, control = TRUE)
  )
}
expected_species <- expected_species %>%
  inner_join(synonyms, by = c("species" = "synonym"))

# read sequence counts summary from FASTQ sampling
counts <- read_csv(counts_file, col_types = "cii")
missing_columns <- setdiff(c("id", "read", "sampled"), colnames(counts))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", counts_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# read bowtie alignment output
alignments <- read_tsv(alignments_file, col_types = "cccciccic")
alignment_columns <- c("genome", "read", "strand", "chr", "pos", "sequence", "quality", "num", "mismatches")
missing_columns <- setdiff(alignment_columns, colnames(alignments))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# calculate sequence lengths, count mismatches and separate sample name from read id
alignments <- alignments %>%
  mutate(length = nchar(sequence)) %>%
  mutate(number_of_mismatches = str_count(mismatches, ",") + 1L) %>%
  mutate(number_of_mismatches = as.integer(replace_na(number_of_mismatches, 0L))) %>%
  separate(read, into = c("id", "read"), sep = ":") %>%
  select(id, read, length, genome, mismatches = number_of_mismatches)

# mark alignments to the expected genome(s)
alignments <- alignments %>%
  left_join(transmute(expected_species, id, genome, expected = TRUE, control), by = c("id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  mutate(control = replace_na(control, FALSE))

# alignment summary
aligned_summary <- alignments %>%
  group_by(id, genome, expected, control) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(counts, id, total = sampled), by = "id") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(id, genome, expected, control, count, percentage, error_rate)

# filter best alignments
if (nrow(alignments) > 0) {
  alignments <- alignments %>%
    group_by(id, read) %>%
    filter(mismatches <= min(mismatches)) %>%
    ungroup()
}

best_aligned_counts <- count(alignments, id, genome, name = "count")

# assign reads to genomes/species
assigned_alignments <- alignments %>%
  left_join(best_aligned_counts, by = c("id", "genome")) %>%
  group_by(id, read) %>%
  arrange(desc(expected), desc(count)) %>%
  slice(1) %>%
  ungroup()

assigned_summary <- assigned_alignments %>%
  group_by(id, genome, expected, control) %>%
  summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
  left_join(select(counts, id, total = sampled), by = "id") %>%
  mutate(percentage = round(100 * count / total, digits = 1)) %>%
  select(id, genome, expected, control, count, percentage, error_rate)

# unmapped read counts
unmapped_counts <- assigned_alignments %>%
  count(id, name = "aligned") %>%
  full_join(counts, by = "id") %>%
  mutate(aligned = replace_na(aligned, 0)) %>%
  mutate(count = sampled - aligned) %>%
  mutate(percentage = round(100 * count / sampled, digits = 1)) %>%
  select(id, count, percentage)

# summary
summary <- aligned_summary %>%
  select(id, genome, expected, control, aligned = count, `aligned %` = percentage, `aligned error rate` = error_rate) %>%
  left_join(select(assigned_summary, id, genome, assigned = count, `assigned %` = percentage, `assigned error rate` = error_rate), by = c("id", "genome")) %>%
  bind_rows(transmute(unmapped_counts, id, genome = "unmapped", aligned = count, `aligned %` = percentage, assigned = count, `assigned %` = percentage)) %>%
  mutate(assigned = replace_na(assigned, 0L)) %>%
  mutate(`assigned %` = replace_na(`assigned %`, 0.0))

summary <- samples %>%
  select(id, species) %>%
  left_join(summary, by = "id") %>%
  mutate(id = as_factor(id)) %>%
  arrange(id, is.na(expected), desc(assigned))

summary %>%
  mutate(expected = ifelse(expected, "yes", "no")) %>%
  mutate(control = ifelse(control, "yes", "no")) %>%
  write_csv(summary_file, na = "")

# stacked bar plot
plot <- summary %>%
  left_join(select(counts, id, `reads` = read, sampled), by = "id") %>%
  mutate(category = ifelse(expected, "expected species", "contaminant")) %>%
  mutate(category = ifelse(control, "control", category)) %>%
  mutate(category = replace_na(category, "unmapped")) %>%
  mutate(category = factor(category, levels = rev(c("expected species", "control", "contaminant", "unmapped")))) %>%
  mutate(id = fct_rev(id)) %>%
  mutate(reads = (assigned / sampled) * reads / 1e6) %>%
  mutate(transparency = pmax(1.5 - `aligned error rate`, 0, na.rm = TRUE)) %>%
  mutate(transparency = ifelse(category == "unmapped", max(transparency), transparency)) %>%
  arrange(id, category, desc(reads)) %>%
  mutate(group = row_number()) %>%
  ggplot(aes(x = id, y = reads, group = group, fill = category, alpha = transparency)) +
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

pdf(plot_file, width = 8, height = 1 + nrow(distinct(summary, id)) * 0.55)
print(plot)
dev.off()
