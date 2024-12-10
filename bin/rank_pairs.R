#!/usr/bin/env Rscript
options(readr.show_col_types=FALSE)
library(yyjsonr) # Fast json parser
library(dplyr)
library(tidyr)
library(readr)


args <- commandArgs(trailingOnly = TRUE)
# args <- c("results/cc_lotp_pairs", "references/CCNA_ref.csv", "references/cc_uniprot_mappings.csv")

result_root <- args[1]
ranked_result_dir <- file.path(result_root, "ranked_results")
if (!dir.exists(ranked_result_dir))
  dir.create(ranked_result_dir, recursive = TRUE)

# Read reference files
ccna_ref <- read_csv(as.character(args[2]))
uniprot_ref <- read_csv(as.character(args[3]))

# List rank_001 json file names with full paths
json_ls <- list.files(
  path = file.path(result_root, "colabfold/webserver"),
  pattern = ".*rank_001.*json",
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
matches <- regmatches(basename(json_ls), regexec(pattern, basename(json_ls)))
# Combine json file name components in a data.frame
matches <- matches |>
  do.call(rbind, args = _) |>
  as.data.frame() |>
  setNames(c(
    "filename",
    "paired_uniprotID",
    "rank",
    "model_type",
    "modelID",
    "seedID"
  ))

matches <- matches |>
  mutate(iptm = iptm_vals) |>
  summarise(iptm = mean(iptm),
            .by = c(paired_uniprotID, rank, model_type)) |> # Averages iptm values if there are more than 1 rank1 json file
  separate_wider_delim(cols = paired_uniprotID,
                       delim = "-",
                       names = c("id1", "id2"),
                       cols_remove = FALSE) |>
  left_join(uniprot_ref, by = join_by(id2 == uniprotID)) |>
  rename(prey_locus_tag = locus_tag) |>
  left_join(uniprot_ref, by = join_by(id1 == uniprotID)) |>
  rename(bait_locus_tag = locus_tag) |>
  select(-id1, -id2) |>
  left_join(ccna_ref, by = join_by(prey_locus_tag == locus_tag)) |>
  arrange(desc(iptm)) |>
  relocate(all_of(c("bait_locus_tag", "prey_locus_tag", "paired_uniprotID", "rank", "model_type", "iptm")))

# Export results to csv
write_csv(matches, file.path(ranked_result_dir, "ranked_results.csv"))
