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

  make_option(c("--output-pdf"), dest = "output_pdf_file",
              help = "Output stacked bar plot PDF file"),

  make_option(c("--output-png"), dest = "output_png_file",
              help = "Output stacked bar plot PNG file"),

  make_option(c("--output-svg"), dest = "output_svg_file",
              help = "Output stacked bar plot SVG file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

summary_file <- opt$summary_file
alignment_summary_file <- opt$alignment_summary_file
output_pdf_file <- opt$output_pdf_file
output_png_file <- opt$output_png_file
output_svg_file <- opt$output_svg_file

if (is.null(summary_file)) stop("Summary file must be specified")
if (is.null(alignment_summary_file)) stop("Alignment summary file must be specified")

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(forcats)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

# genome alignment summaries
alignment_summary <- read_csv(
  alignment_summary_file,
  col_types = cols_only(
    id = col_character(),
    genome = col_character(),
    expected = col_character(),
    control = col_character(),
    assigned = col_double(),
    `assigned error rate` = col_double()
  )
)

alignment_summary <- alignment_summary %>%
  rename(user_id = id) %>%
  mutate(expected = expected == "yes") %>%
  mutate(control = control == "yes")

sample_ids <- alignment_summary %>%
  distinct(user_id) %>%
  mutate(id = row_number())

alignment_summary <- alignment_summary %>%
  left_join(sample_ids, by = "user_id") %>%
  select(id, everything())

# data set summaries
summary <- read_csv(
  summary_file,
  col_types = cols_only(
    id = col_character(),
    sequences = col_double(),
    sampled = col_double(),
    adapter = col_double()
  )
)

summary <- summary %>%
  rename(user_id = id) %>%
  left_join(sample_ids, by = "user_id") %>%
  select(id, everything())

counts <- select(summary, id, sequences, sampled)

lower_error_rate <- 0.25
upper_error_rate <- 1.5

plot_data <- alignment_summary %>%
  select(id, user_id, reference = genome, count = assigned, error_rate = `assigned error rate`, expected, control) %>%
  filter(count > 0) %>%
  mutate(category = ifelse(expected, "expected species", "unexpected/contaminant")) %>%
  mutate(category = ifelse(control, "control", category)) %>%
  mutate(category = replace_na(as.character(category), "unmapped")) %>%  # category will be logical if there are no rows with non-zero counts
  select(-expected, -control) %>%
  mutate(type = "genome") %>%
  bind_rows(transmute(summary, id, user_id, reference = "adapter", count = adapter, error_rate = lower_error_rate, category = "adapter", type = "adapter")) %>%
  left_join(counts, by = "id") %>%
  mutate(scaled_count = ifelse(sampled == 0, 0, (count / sampled) * sequences / 1e6)) %>%
  select(-sampled, -sequences) %>%
  arrange(id) %>%
  mutate(user_id = as_factor(user_id)) %>%
  mutate(type = factor(type, levels = c("adapter", "genome"))) %>%
  mutate(category = factor(category, levels = rev(c("expected species", "control", "unexpected/contaminant", "unmapped", "adapter")))) %>%
  mutate(bounded_error_rate = pmin(error_rate, upper_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = pmax(bounded_error_rate, lower_error_rate, na.rm = TRUE)) %>%
  mutate(bounded_error_rate = ifelse(category == "unmapped", lower_error_rate, bounded_error_rate)) %>%
  arrange(user_id, type, category, count) %>%
  mutate(group = row_number()) %>%
  mutate(bar_width = ifelse(type == "adapter", 0.5, 1))

plot <- ggplot(plot_data, aes(x = type, y = scaled_count, group = group, fill = category, alpha = -bounded_error_rate)) +
  geom_col(colour = "grey30", width = plot_data$bar_width) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = 0), labels = scales::comma) +
  scale_fill_manual(values = c("magenta", "grey97", "red", "goldenrod1", "green"), drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_alpha(guide = "none") +
  labs(y = "sequence reads (millions)") +
  coord_flip() +
  facet_grid(rows = vars(user_id), switch = "y") +
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
height <- 0.7 * (1.5 + nrow(sample_ids))

if (!is.null(output_pdf_file)) {
  ggsave(output_pdf_file, plot, width = width, height = height, limitsize = FALSE)
}

if (!is.null(output_png_file)) {
  ggsave(output_png_file, plot, width = width, height = height, limitsize = FALSE)
}

if (!is.null(output_svg_file)) {
  ggsave(output_svg_file, plot, width = width, height = height, limitsize = FALSE)
}
