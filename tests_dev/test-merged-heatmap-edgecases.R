# Reproduce signed cancellation, negative-only and zero merged heatmap cases.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
source("tests_dev/test-merge-spatialcellchat.R", local = FALSE)
out.dir <- "tests_dev/figures/merged/heatmap-edgecases"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
probe <- function(case) {
  obj <- merged
  first <- obj@misc$.datasets$sample1$net$group$count
  second <- obj@misc$.datasets$sample2$net$group$count
  first[,] <- 0; second[,] <- 0
  if (case == "cancel") { first[1, 2] <- 1; second[1, 1] <- 1 }
  if (case == "balanced_isolates") {
    first[1, 2] <- 1; first[2, 1] <- 1
    second[1, 1] <- 1; second[2, 2] <- 1
  }
  if (case == "negative") first[1, 2] <- 1
  obj@misc$.datasets$sample1$net$group$count <- first
  obj@misc$.datasets$sample2$net$group$count <- second
  delta <- as.matrix(second - first)
  file <- file.path(out.dir, paste0(case, ".pdf"))
  message <- tryCatch({
    grDevices::pdf(file, width = 8, height = 6)
    result <- tryCatch(netVisual_heatmap(obj, comparison = c(1, 2), measure = "count",
                                        remove.isolate = case %in% c("balanced_isolates", "zero_isolates")),
                       finally = grDevices::dev.off())
    plotted <- result@ht_list[[1L]]@matrix
    if (case == "balanced_isolates") {
      stopifnot(identical(rownames(plotted), rownames(delta)[1:2]),
                identical(colnames(plotted), colnames(delta)[1:2]))
      delta <- delta[1:2, 1:2, drop = FALSE]
    }
    stopifnot(all(plotted[delta != 0] == delta[delta != 0]))
    ""
  }, error = conditionMessage)
  if (nzchar(message) && file.exists(file)) unlink(file)
  data.frame(case = case, sum = sum(delta), nonzero = sum(delta != 0),
             status = if (nzchar(message)) "ERROR" else "PDF-WRITTEN", detail = message,
             stringsAsFactors = FALSE)
}
summary <- do.call(rbind, lapply(c("cancel", "negative", "balanced_isolates", "zero", "zero_isolates"), probe))
print(summary)
write.csv(summary, file.path(out.dir, "summary.csv"), row.names = FALSE)
if (summary$status[1] != "PDF-WRITTEN" || summary$status[2] != "PDF-WRITTEN" ||
    summary$status[3] != "PDF-WRITTEN" || any(summary$status[4:5] != "ERROR") ||
    !all(grepl("No significant cell-cell communication", summary$detail[4:5], fixed = TRUE)))
  quit(status = 1L)
