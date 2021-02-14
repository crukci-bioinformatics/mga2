args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
{
  stop('Usage: Rscript check_inputs.R samples_file genomes_file output_file')
}

samples_file <- args[1]
genomes_file <- args[2]
output_file <- args[3]

suppressPackageStartupMessages(library(tidyverse))

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
missing_columns <- setdiff(c("id", "fastq", "species", "control"), colnames(samples))
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
  mutate(synonym = str_trim(synonym)) %>%
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

samples %>%
  left_join(synonyms, by = c("species" = "synonym")) %>%
  left_join(select(synonyms, control = synonym, control_genome = genome), by = "control") %>%
  mutate(id_prefix = row_number()) %>%
  write_csv(output_file)
