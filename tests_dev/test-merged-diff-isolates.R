# Check signed isolate pruning in actual merged differential network rendering.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
source("tests_dev/test-merge-spatialcellchat.R", local = FALSE)
out.dir <- "tests_dev/figures/merged/diff-isolate-edgecases"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
original.graph <- graph_from_adjacency_matrix
captured.net <- NULL
graph_from_adjacency_matrix <- function(adjmatrix, ...) {
  captured.net <<- adjmatrix
  original.graph(adjmatrix, ...)
}
probe <- function(case) {
  obj <- merged
  first <- obj@misc$.datasets$sample1$net$group$count
  second <- obj@misc$.datasets$sample2$net$group$count
  first[,] <- 0; second[,] <- 0
  if (case == "balanced") {
    first[1, 2] <- 1; first[2, 1] <- 1
    second[1, 1] <- 1; second[2, 2] <- 1
    retained <- rownames(first)[1:2]
  } else if (case == "no_isolates") {
    diag(second) <- 1
    retained <- rownames(first)
  } else if (case == "one_node") {
    second[1, 1] <- 1
    retained <- rownames(first)[1]
  } else if (case == "empty") {
    retained <- character()
  } else stop("unknown case")
  obj@misc$.datasets$sample1$net$group$count <- first
  obj@misc$.datasets$sample2$net$group$count <- second
  expected <- as.matrix(second - first)[retained, retained, drop = FALSE]
  file <- file.path(out.dir, paste0(case, ".pdf"))
  captured.net <<- NULL
  message <- tryCatch({
    grDevices::pdf(file, width = 8, height = 6)
    tryCatch(netVisual_diffInteraction(obj, comparison = c(1, 2), measure = "count", remove.isolate = TRUE),
             finally = grDevices::dev.off())
    stopifnot(identical(as.matrix(captured.net), expected))
    ""
  }, error = conditionMessage)
  if (nzchar(message) && file.exists(file)) unlink(file)
  data.frame(case = case, retained = length(retained), status = if (nzchar(message)) "ERROR" else "PDF-WRITTEN", detail = message)
}
summary <- do.call(rbind, lapply(c("balanced", "no_isolates", "one_node", "empty"), probe))
write.csv(summary, file.path(out.dir, "summary.csv"), row.names = FALSE)
print(summary)
if (any(summary$status[1:3] != "PDF-WRITTEN") || summary$status[4] != "ERROR" ||
    !grepl("No nonzero differential interactions", summary$detail[4], fixed = TRUE)) quit(status = 1L)
