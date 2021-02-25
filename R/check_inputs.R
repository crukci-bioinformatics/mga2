args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
{
  stop('Usage: Rscript check_inputs.R samples_file genome_details_file bowtie_index_list_file output_samples_file output_genomes_file')
}

samples_file <- args[1]
genome_details_file <- args[2]
bowtie_index_list_file <- args[3]
output_samples_file <- args[4]
output_genomes_file <- args[5]

suppressPackageStartupMessages(library(tidyverse))

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))

expected_columns <- c("id", "fastq", "species", "controls")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

if (nrow(samples) == 0) {
  stop("empty sample sheet: ", samples_file)
}

if (nrow(filter(samples, is.na(samples$id))) > 0) {
  stop("missing ids found in ", samples_file)
}

duplicates <- samples %>%
  count(id) %>%
  filter(n > 1) %>%
  pull(id)
if (length(duplicates) > 0) {
  stop("duplicate ids found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# read list of genomes with bowtie indexes
bowtie_indexes <- read_tsv(bowtie_index_list_file, col_names = "id", col_types = "c") %>%
  mutate(id_lower_case = str_to_lower(id))

invalid_names <- bowtie_indexes %>%
  filter(str_detect(id, "\\|")) %>%
  arrange(id) %>%
  pull(id)
if (length(invalid_names) > 0) {
  stop("the following bowtie index names contain the '|' character: '", str_c(invalid_names, collapse = "', '"), "'")
}

duplicates <- bowtie_indexes %>%
  add_count(id_lower_case) %>%
  filter(n > 1) %>%
  arrange(id_lower_case, id) %>%
  pull(id)
if (length(duplicates) > 0) {
  stop("the following bowtie index names differ only by case: '", str_c(duplicates, collapse = "', '"), "'")
}

# read genome details file
genome_details <- read_csv(genome_details_file, col_types = cols(.default = col_character()))

expected_columns <- c("genome", "species", "synonyms")
missing_columns <- setdiff(expected_columns, colnames(genome_details))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", genome_details_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

invalid_names <- genome_details %>%
  filter(str_detect(genome, "\\|")) %>%
  arrange(genome) %>%
  pull(genome)
if (length(invalid_names) > 0) {
  stop("the following genome names contain the '|' character: '", str_c(invalid_names, collapse = "', '"), "'")
}

duplicates <- genome_details %>%
  count(genome) %>%
  filter(n > 1) %>%
  pull(genome)
if (length(duplicates) > 0) {
  stop("duplicate genome ids found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

genome_details <- genome_details %>%
  select(all_of(expected_columns)) %>%
  mutate(genome_lower_case = str_to_lower(genome))

if (nrow(filter(genome_details, is.na(genome_details$genome))) > 0) {
  stop("missing genome ids found in ", genome_details_file)
}

duplicates <- genome_details %>%
  add_count(genome_lower_case) %>%
  filter(n > 1) %>%
  arrange(genome_lower_case, genome) %>%
  pull(genome)
if (length(duplicates) > 0) {
  stop("genome ids differing only by case found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# match genome details to bowtie indexes
# note that the bowtie prefix takes precedence over the genome id
# in the genome details file if these differ by case
genome_details <- genome_details %>%
  full_join(bowtie_indexes, by = c("genome_lower_case" = "id_lower_case")) %>%
  select(-genome_lower_case) %>%
  mutate(has_index = !is.na(id)) %>%
  mutate(genome = ifelse(has_index, id, genome)) %>%
  select(-id)

# create synonym lookup table
synonyms <- bind_rows(
  genome_details %>%
    select(genome = genome, synonym = genome),
  genome_details %>%
    select(genome, synonym = species),
  genome_details %>%
    select(genome, synonym = synonyms) %>%
    separate_rows(synonym, sep = "\\|")
) %>%
  filter(!is.na(synonym)) %>%
  mutate(synonym = str_trim(synonym)) %>%
  filter(nchar(synonym) > 0) %>%
  mutate(synonym = str_to_lower(synonym)) %>%
  distinct()

duplicates <- synonyms %>%
  count(synonym, name = "count") %>%
  filter(count > 1) %>%
  pull(synonym)

if (length(duplicates) > 0) {
  stop("synonyms matching (case-insensitive) more than one genome found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# combine genome details and synonyms tables
synonyms <- genome_details %>%
  select(-synonyms) %>%
  left_join(synonyms, by = "genome")

# add unique numeric id to be used in labeling sampled reads
samples <- samples %>%
  mutate(rownum = row_number()) %>%
  select(rownum, everything())

# find matching genomes
species <- samples %>%
  select(rownum, species) %>%
  separate_rows(species, sep = "\\|") %>%
  drop_na() %>%
  mutate(synonym = str_to_lower(str_trim(species)))

non_matching <- species %>%
  anti_join(synonyms, by = "synonym") %>%
  distinct(species) %>%
  pull(species)
if (length(non_matching) > 0) {
  message("No genome matches found for the following species: '", str_c(non_matching, collapse = "', '"), "'")
}

genomes <- species %>%
  inner_join(synonyms, by = "synonym") %>%
  select(rownum, genome)

# find matching control genomes
controls <- samples %>%
  select(rownum, controls) %>%
  separate_rows(controls, sep = "\\|") %>%
  drop_na() %>%
  mutate(synonym = str_to_lower(str_trim(controls)))

non_matching <- controls %>%
  anti_join(synonyms, by = "synonym") %>%
  distinct(controls) %>%
  pull(controls)
if (length(non_matching) > 0) {
  message("No genome matches found for the following controls: '", str_c(non_matching, collapse = "', '"), "'")
}

control_genomes <- controls %>%
  inner_join(synonyms, by = "synonym") %>%
  select(rownum, genome)

# remove genome matches if these also match as a control
genomes <- anti_join(genomes, control_genomes, by = c("rownum", "genome"))

# collapse multiple genomes for the same sample/dataset
genomes <- genomes %>%
  group_by(rownum) %>%
  summarize(genomes = str_c(genome, collapse = "|"), .groups = "drop")

control_genomes <- control_genomes %>%
  group_by(rownum) %>%
  summarize(control_genomes = str_c(genome, collapse = "|"), .groups = "drop")

# append matching genomes and control genomes
samples <- samples %>%
  left_join(genomes, by = "rownum") %>%
  left_join(control_genomes, by = "rownum")

# write new sample sheet
write_csv(samples, output_samples_file, na = "")

# write genomes output file that includes matching species for each genome
genome_details %>%
  filter(has_index) %>%
  select(genome, species) %>%
  write_csv(output_genomes_file, na = "")


