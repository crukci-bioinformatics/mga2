args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
{
  stop('Usage: Rscript check_inputs.R samples_file genomes_file genome_list_file output_file')
}

args <- c("samplesheet.csv", "genomes.csv", "genome_list.txt", "samples.validated.csv")

samples_file <- args[1]
genomes_file <- args[2]
genome_list_file <- args[3]
output_file <- args[4]

suppressPackageStartupMessages(library(tidyverse))

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
missing_columns <- setdiff(c("id", "fastq", "species", "control"), colnames(samples))
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

# read genomes file
genomes <- read_csv(genomes_file, col_types = cols(.default = col_character()))
missing_columns <- setdiff(c("genome", "species", "synonyms"), colnames(genomes))
if (length(missing_columns) > 0) {
  stop("Missing columns in ", genomes_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

# read list of genomes with bowtie indexes
genome_list <- read_tsv(genome_list_file, col_names = "genome", col_types = "c")

# add missing genomes in the genomes file for which we have bowtie indexes
genomes <- full_join(genomes, genome_list, by = "genome")

# create synonym lookup table
synonyms <- bind_rows(
  genomes %>%
    select(genome = genome, synonym = genome),
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

# match species and controls to genomes
samples <- samples %>%
  left_join(synonyms, by = c("species" = "synonym")) %>%
  left_join(select(synonyms, control = synonym, `control genome` = genome), by = "control")

# check for species that don't match to a genome
non_matching_species <- samples %>%
  filter(!is.na(species)) %>%
  filter(is.na(genome)) %>%
  distinct(species) %>%
  pull(species)
if (length(non_matching_species) > 0) {
  stop("no genome match for the following species: '", str_c(non_matching_species, collapse = "', '"), "'")
}

# check for controls that don't match to a genome
non_matching_control <- samples %>%
  filter(!is.na(control)) %>%
  filter(is.na(`control genome`)) %>%
  distinct(control) %>%
  pull(control)
if (length(non_matching_control) > 0) {
  stop("no genome match for the following controls: '", str_c(non_matching_control, collapse = "', '"), "'")
}

# add unique numeric identifier to be used in labeling sampled reads
samples <- mutate(samples, id_prefix = row_number())

# write new sample sheet
write_csv(samples, output_file)
