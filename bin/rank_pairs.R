#!/usr/bin/env Rscript
options(readr.show_col_types=FALSE)
library(yyjsonr) # Fast json parser
library(dplyr)
library(tidyr)
library(readr)


args <- commandArgs(trailingOnly = TRUE)
# args <- c("tests/uniprotkb_proteome_UP000001364_cc.tsv")

# Read reference files
ref_file <- read_tsv(as.character(args[1]))

# List rank_001 json file names with full paths
json_ls <- list.files(
  path = ".",
  pattern = ".*scores_rank_001.*json",
  full.names = TRUE
)
# Collect iptm values from all json files in a numeric vector
iptm_vals <- vapply(json_ls,
                    \(x) read_json_file(x)[["iptm"]],
                    FUN.VALUE = double(1),
                    USE.NAMES = FALSE)

# Regular expression to capture the components of the json files names
pattern <- "^([A-Z0-9-]+)_scores_(rank_\\d{3})_(.*)_(model_\\d{1})_(seed_\\d{3})\\.json$"
# Apply the regex to extract the components
matches <- regmatches(basename(json_ls), regexec(pattern, basename(json_ls))) |>
  do.call(rbind, args = _) |>
  as.data.frame() |>
  setNames(c(
    "filename",
    "paired_uniprotID",
    "rank",
    "model_type",
    "modelID",
    "seedID"
  )) |>
  mutate(iptm = iptm_vals) |>
  summarise(iptm = mean(iptm),
            .by = c(paired_uniprotID, rank, model_type)) |> # Averages iptm values if there are more than 1 rank1 json file
  separate_wider_delim(cols = paired_uniprotID,
                       delim = "-",
                       names = c("id1", "id2"),
                       cols_remove = FALSE) |>
  left_join(ref_file, by = join_by(id2 == Entry)) |>
  left_join(ref_file, by = join_by(id1 == Entry),
            suffix = c("_prey", "_bait")) |>
  select(-id1, -id2) |>
  arrange(desc(iptm))

# Export results to csv
write_tsv(matches, "ranked_results.tsv")
