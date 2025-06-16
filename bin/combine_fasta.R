#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(
    c("--acc_file"),
    type = "character",
    default = NULL,
    help = "Path to tsv file with accession numbers.\n Headers: geneID, bait"
  ),
  make_option(
    c("--proj_dir"),
    type = "character",
    default = NULL,
    help = "Path to project directory"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Pull rest of the libraries
library(seqinr)
library(httr2)

get_fasta <- function(base_url, header, bait_status) {
  sp <- request(base_url) |>
    req_perform() |>
    resp_body_string()

  # Split into lines
  lines <- unlist(strsplit(sp, "\n"))
  # Get sequence
  sequence <- paste(lines[-1], collapse = "")

  # Create dataframe
  data.frame(
    fasta_names = header,
    fasta_sequence = sequence,
    bait_status = bait_status,
    stringsAsFactors = FALSE
  )
}

acclist <- read.table(
  opt$acc_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
# Identify fasta files and uniprot accessions and bait status
fasta_files <- grepl("\\.(fasta|fa)$", acclist$Entry, ignore.case = TRUE)
acc_fasta <- acclist$Entry[fasta_files]
acc_fasta_bait <- acclist$bait[fasta_files]
acc_uniprot <- acclist$Entry[!fasta_files]
acc_uniprot_bait <- acclist$bait[!fasta_files]

# Check if paths are absolute or relative
is_absolute_path <- function(path) {
  # On Unix-like systems, absolute paths start with "/"
  # On Windows, they start with a drive letter (e.g., "C:/") or "\\"
  grepl("^(/|[A-Za-z]:/|\\\\)", path)
}

# Read fasta files into a dataframe with bait status
fasta_df <- data.frame()
if (length(acc_fasta) > 0) {
  # Copy files to workdir
  for (f in acc_fasta) {
    if (is_absolute_path(f)) {
      file.copy(f, ".")
  } else {
      file.copy(file.path(opt$proj_dir, f), ".")
  }
  }

  fasta_df <- do.call(
    rbind,
    lapply(seq_along(acc_fasta), function(i) {
      fasta_content <- read.fasta(
        basename(acc_fasta[i]),
        seqtype = "AA",
        as.string = TRUE,
        set.attributes = FALSE
      )
      data.frame(
        fasta_names = names(fasta_content),
        fasta_sequence = unlist(fasta_content, use.names = FALSE),
        bait_status = acc_fasta_bait[i],
        stringsAsFactors = FALSE
      )
    })
  )
  # Clean the input fastas from work dir
  file.remove(basename(acc_fasta))
}

# For non-fasta entries, use get_fasta to fetch from Uniprot
if (length(acc_uniprot) > 0) {
  uniprot_df <- do.call(
    rbind,
    lapply(seq_along(acc_uniprot), function(i) {
      get_fasta(
        sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc_uniprot[i]),
        acc_uniprot[i],
        acc_uniprot_bait[i]
      )
    })
  )
}

fasta_df <- rbind(fasta_df, uniprot_df)

# Combine bait-prey pairs
combs <- expand.grid(
  fasta_df[fasta_df$bait_status == TRUE, "fasta_names"],
  fasta_df[fasta_df$bait_status == FALSE, "fasta_names"],
  stringsAsFactors = FALSE
)
combs$combined <- paste(combs$Var1, combs$Var2, sep = "-")

# Function to pull sequences from the dataframe based on fasta_names
pull_sequence <- function(df, fasta_name) {
  unique(df$fasta_sequence[df$fasta_names == fasta_name])
}

fasta_out <- mapply(
  \(x, y) {
    paste(pull_sequence(fasta_df, x), pull_sequence(fasta_df, y), sep = ":")
  },
  combs$Var1,
  combs$Var2,
  SIMPLIFY = FALSE
) |>
  setNames(combs$combined)

# Write pair results to fasta files
for (i in seq_along(fasta_out)) {
  write.fasta(
    fasta_out[[i]],
    names = names(fasta_out)[i],
    file.out = paste0(names(fasta_out)[i], ".fasta"),
    as.string = TRUE
  )
}
