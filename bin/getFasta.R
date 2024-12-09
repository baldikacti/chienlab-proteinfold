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
