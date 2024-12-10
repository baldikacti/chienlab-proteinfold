#!/usr/bin/env Rscript
library(optparse)
library(seqinr)

option_list <- list(
  make_option(c("--acc_file"),
              type = "character", default = NULL,
              help = "Path to csv file with accession numbers.\n Headers: geneID, bait"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

getFasta <- function(base_url, outfile) {
  library(httr2)
  pb <- progress::progress_bar$new(
    format = "Downloading Fasta [:bar] :percent",
    total = length(base_url),
    clear = FALSE, width = 60
  )
  for (i in seq_along(base_url)) {
    tryCatch(
      {
        request(base_url[i]) |>
          req_perform() |>
          resp_body_string() |>
          write.table(
            file = paste0(outfile, ".fasta"),
            quote = F, row.names = F, col.names = F, append = T
          )
        message(paste0("Downloaded: ", basename(base_url), "\n"))
      },
      error = function(e) {
        message(paste("Encountered error processing:", base_url[i]))
        message(conditionMessage(e))
        NA
      }
    )

    Sys.sleep(1)
    pb$tick()
  }
}

acclist <- read.table(opt$acc_file, header = T, sep = "\t")
query_url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acclist[[1]])

getFasta(query_url, outfile = "query")

# Reads the fasta file and creates combinations of proteins
fasta <- read.fasta("query.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE) |>
  setNames(acclist[[1]])
combs <- expand.grid(acclist[acclist[[2]] == TRUE, 1], acclist[acclist[[2]] == FALSE, 1], stringsAsFactors = FALSE)
combs$combined <- paste(combs$Var1, combs$Var2, sep = "-")
fast_out <- mapply(\(x,y) {paste(fasta[[x]], fasta[[y]], sep = ":")}, combs$Var1, combs$Var2, SIMPLIFY = FALSE) |>
  setNames(combs$combined)

for (i in seq_along(fast_out)) {
  write.fasta(fast_out[[i]], names = names(fast_out)[i], file.out = paste0(names(fast_out)[i], ".fasta"), as.string = TRUE)
}

#samplesheet <- data.frame(
#  sequence = names(fast_out),
#  fasta = paste0(names(fast_out), ".fasta")
#)

#write.csv(samplesheet, opt$samplesheet, row.names = FALSE, quote = FALSE)

# Clean up
file.remove("query.fasta")

