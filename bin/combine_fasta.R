#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--acc_file"),
              type = "character", default = "tests/acclist.csv",
              help = "Path to csv file with accession numbers.\n Headers: geneID, bait"
  ),
  make_option(c("--fasta_dir"),
              type = "character", default = "./fasta",
              help = "Directory to export fasta files"
  ),
  make_option(c("--samplesheet"),
              type = "character", default = "tests/samplesheet.csv",
              help = "CSV file to export processed accession file."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

library(seqinr)
source("bin/getFasta.R")

# Initialize directories and clean preivous query.fasta
if(!dir.exists(opt$fasta_dir)) dir.create(opt$fasta_dir, recursive = TRUE)
if(!dir.exists(dirname(opt$samplesheet))) dir.create(dirname(opt$samplesheet), recursive = TRUE)
if(file.exists(file.path(opt$fasta_dir, "query.fasta"))) file.remove(file.path(opt$fasta_dir, "query.fasta"))

acclist <- read.csv(opt$acc_file, header = T)
query_url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acclist[[1]])

getFasta(query_url, outfile = file.path(opt$fasta_dir, "query"))

# Reads the fasta file and creates combinations of proteins
fasta <- read.fasta(file.path(opt$fasta_dir, "query.fasta"), seqtype = "AA", as.string = TRUE, set.attributes = FALSE) |>
  setNames(acclist[[1]])
combs <- expand.grid(acclist[acclist[[2]] == TRUE, "geneID"], acclist[acclist[[2]] == FALSE, 1], stringsAsFactors = FALSE)
combs$combined <- paste(combs$Var1, combs$Var2, sep = "-")
fast_out <- mapply(\(x,y) {paste(fasta[[x]], fasta[[y]], sep = ":")}, combs$Var1, combs$Var2, SIMPLIFY = FALSE) |>
  setNames(combs$combined)

for (i in seq_along(fast_out)) {
  write.fasta(fast_out[[i]], names = names(fast_out)[i], file.out = file.path(opt$fasta_dir, paste0(names(fast_out)[i], ".fasta")), as.string = TRUE)
}

samplesheet <- data.frame(
  sequence = names(fast_out),
  fasta = file.path(opt$fasta_dir, paste0(names(fast_out), ".fasta"))
)

write.csv(samplesheet, opt$samplesheet, row.names = FALSE, quote = FALSE)

# Clean up
file.remove(file.path(opt$fasta_dir, "query.fasta"))

