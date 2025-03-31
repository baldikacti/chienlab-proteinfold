#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--acc_file"),
    type = "character", default = NULL,
    help = "Path to tsv file with accession numbers.\n Headers: geneID, bait"
  ),
  make_option(c("--proj_dir"),
    type = "character", default = NULL,
    help = "Path to project directory"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

library(seqinr)
getFasta <- function(base_url, outfile) {
  library(httr2)
  pb <- progress::progress_bar$new(
    format = "Downloading Fasta [:bar] :percent",
    total = length(base_url),
    clear = FALSE, width = 60
  )
  for (url in base_url) {
    tryCatch(
      {
        request(url) |>
          req_perform() |>
          resp_body_string() |>
          write.table(
            file = paste0(outfile, ".fasta"),
            quote = F, row.names = F, col.names = F, append = T
          )
        message(paste0("Downloaded: ", basename(url), "\n"))
      },
      error = function(e) {
        message(paste("Encountered error processing:", base_url[i]))
        message(conditionMessage(e))
        NA
      }
    )

    Sys.sleep(0.1) # To not spam uniprot
    pb$tick()
  }
}

acclist <- read.table(opt$acc_file, header = T, sep = "\t")
fasta_files <- grepl("\\.fasta$", acclist[[1]])
acc_fasta <- acclist[fasta_files, 1]
acc_fasta_names <- sub(".fasta", "", basename(acc_fasta))

acc_uniprot <- acclist[!fasta_files, 1]

query_url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc_uniprot)

getFasta(query_url, outfile = "query")

# Reads the fasta files and creates combinations of proteins
fasta <- read.fasta("query.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE) |>
  setNames(acc_uniprot)
if (sum(fasta_files) > 0) {
  # Copy files to workdir
  file.copy(file.path(opt$proj_dir, acc_fasta), ".")
  fasta_file <- lapply(basename(acc_fasta), \(x) read.fasta(x, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)) |>
    unlist(recursive = FALSE)
  fasta <- c(fasta, fasta_file)
  # Sanitize input names if it contains file paths
  acclist[fasta_files, 1] <- acc_fasta_names
  file.remove(basename(acc_fasta))
}

# Combine bait-prey pairs
combs <- expand.grid(acclist[acclist[[2]] == TRUE, 1], acclist[acclist[[2]] == FALSE, 1], stringsAsFactors = FALSE)
combs$combined <- paste(combs$Var1, combs$Var2, sep = "-")
fast_out <- mapply(\(x, y) {
  paste(fasta[[x]], fasta[[y]], sep = ":")
}, combs$Var1, combs$Var2, SIMPLIFY = FALSE) |>
  setNames(combs$combined)

# Write pair results to fasta files
for (i in seq_along(fast_out)) {
  write.fasta(fast_out[[i]], names = names(fast_out)[i], file.out = paste0(names(fast_out)[i], ".fasta"), as.string = TRUE)
}

# Clean up
file.remove("query.fasta")
