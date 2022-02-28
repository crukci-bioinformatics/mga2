#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Create stacked bar chart plot.

options(bitmapType='cairo')

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--summary"), dest = "summary_file",
              help = "CSV file containing alignment/assignment summary for each data set"),

  make_option(c("--alignment-summary"), dest = "alignment_summary_file",
              help = "CSV file containing alignment/assignment summary for each genome and data set"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output plot files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

summary_file <- opt$summary_file
alignment_summary_file <- opt$alignment_summary_file
output_prefix <- opt$output_prefix

if (is.null(summary_file)) stop("Summary file must be specified")
if (is.null(alignment_summary_file)) stop("Alignment summary file must be specified")

if (is.null(output_prefix)) output_prefix <- ""

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(forcats)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

output_plot_file_prefix <- str_c(output_prefix, "mga_alignment_summary")

# genome alignment summaries
alignment_summary <- read_csv(alignment_summary_file, show_col_types = FALSE)

required_columns <- c("id", "genome", "expected", "control", "assigned", "assigned error rate")
missing_columns <- setdiff(required_columns, colnames(alignment_summary))
if (length(missing_columns) > 0) {
  stop("missing columns in ", alignment_summary_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

alignment_summary <- alignment_summary %>%
  select(one_of(required_columns)) %>%
  rename(name = id) %>%
  mutate(expected = expected == "yes") %>%
  mutate(control = control == "yes")

sample_names <- alignment_summary %>%
  distinct(name) %>%
  mutate(id = row_number())

alignment_summary <- alignment_summary %>%
  left_join(sample_names, by = "name") %>%
  select(id, everything())

# data set summaries
summary <- read_csv(summary_file, show_col_types = FALSE)

required_columns <- c("id", "sequences", "sampled", "adapter")
missing_columns <- setdiff(required_columns, colnames(summary))
if (length(missing_columns) > 0) {
  stop("missing columns in ", summary_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

summary <- summary %>%
  select(one_of(required_columns)) %>%
  rename(name = id) %>%
  left_join(sample_names, by = "name") %>%
  select(id, everything())

counts <- select(summary, id, sequences, sampled)

lower_error_rate <- 0.25
upper_error_rate <- 1.5

plot_data <- alignment_summary %>%
  select(id, name, reference = genome, count = assigned, error_rate = `assigned error rate`, expected, control) %>%
  filter(count > 0) %>%
  mutate(category = ifelse(expected, "expected species", "unexpected/contaminant")) %>%
  mutate(category = ifelse(control, "control", category)) %>%
  mutate(category = replace_na(category, "unmapped")) %>%
  select(-expected, -control) %>%
  mutate(type = "genome") %>%
  bind_rows(transmute(summary, id, name, reference = "adapter", count = adapter, error_rate = lower_error_rate, category = "adapter", type = "adapter")) %>%
  left_join(counts, by = "id") %>%
  mutate(scaled_count = ifelse(sampled == 0, 0, (count / sampled) * sequences / 1e6)) %>%
  select(-sampled, -sequences) %>%
  arrange(id) %>%
  mutate(name = as_factor(name)) %>%
  mutate(type = factor(type, levels = c("adapter", "genome"))) %>%
  mutate(category = factor(category, levels = rev(c("expected species", "control", "unexpected/contaminant", "unmapped", "adapter")))) %>%
  mutate(bounded_error_rate = pmin(error_rate, upper_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = pmax(bounded_error_rate, lower_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = ifelse(category == "unmapped", lower_error_rate, bounded_error_rate)) %>%
  arrange(name, type, category, count) %>%
  mutate(group = row_number()) %>%
  mutate(bar_width = ifelse(type == "adapter", 0.5, 1))

plot <- ggplot(plot_data, aes(x = type, y = scaled_count, group = group, fill = category, alpha = -bounded_error_rate)) +
  geom_col(colour = "grey30", width = plot_data$bar_width) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = 0), labels = scales::comma) +
  scale_fill_manual(values = c("magenta", "grey97", "red", "goldenrod1", "green"), drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_alpha(guide = "none") +
  labs(y = "sequence reads (millions)") +
  coord_flip() +
  facet_grid(rows = vars(name), switch = "y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.height = unit(0.8, "line"),
    legend.key.width = unit(1.5, "line"),
    legend.text = element_text(margin = margin(r = 1, unit = "line")),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(0.35, "line"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, vjust = 0.8)
  )

width <- 10
height <- 0.7 * (1.5 + nrow(sample_names))

ggsave(str_c(output_plot_file_prefix, ".png"), plot, width = width, height = height)
ggsave(str_c(output_plot_file_prefix, ".pdf"), plot, width = width, height = height)
ggsave(str_c(output_plot_file_prefix, ".svg"), plot, width = width, height = height)

