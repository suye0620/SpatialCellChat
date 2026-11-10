# Render comparison consumers from the same real merge producer used in the regression.
setwd(local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value=TRUE)); if (length(f)) dirname(dirname(normalizePath(f))) else getwd() }))
source("tests_dev/test-merge-spatialcellchat.R", local = FALSE)
out.dir <- "tests_dev/figures/merged/producer-comparison"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
# The two inputs are the simulated fixture and a deliberately modified copy;
# these are examples of producer/consumer behavior, not biological comparisons.
cases <- list(
  compareInteractions = function() compareInteractions(merged, measure = "count"),
  rankNet = function() rankNet(merged, slot.name = "netP", mode = "comparison", comparison = c(1,2), measure = "count", return.data = TRUE)$gg.obj,
  netVisual_barplot = function() netVisual_barplot(merged, comparison = c(1,2), measure = "count"),
  netVisual_diffInteraction = function() netVisual_diffInteraction(merged, comparison = c(1,2), measure = "count"),
  netVisual_heatmap = function() netVisual_heatmap(merged, comparison = c(1,2), measure = "count"),
  netVisual_bubble = function() netVisual_bubble(merged, comparison = c(1,2), signaling = "P1", return.data = TRUE)$gg.obj,
  netAnalysis_diff_signalingRole_scatter = function() netAnalysis_diff_signalingRole_scatter(merged, comparison = c(1,2), do.label = FALSE, show.legend = FALSE),
  netAnalysis_signalingChanges_scatter = function() netAnalysis_signalingChanges_scatter(merged, idents.use = "A", comparison = c(1,2), do.label = FALSE, show.legend = FALSE),
  netVisual_embeddingPairwise = function() netVisual_embeddingPairwise(similarity.merged, slot.name = "netP", comparison = c(1,2), do.label = FALSE, show.legend = FALSE),
  netVisual_embeddingPairwiseZoomIn = function() netVisual_embeddingPairwiseZoomIn(similarity.merged, slot.name = "netP", comparison = c(1,2), pathway.remove = character(), do.label = FALSE, show.legend = FALSE)
)
results <- lapply(names(cases), function(name) {
  out <- file.path(out.dir, paste0(name, ".pdf"))
  detail <- tryCatch({
    grDevices::pdf(out, width = 8, height = 6)
    tryCatch({
      value <- cases[[name]]()
      if (inherits(value, "ggplot")) print(value)
      else if (!name %in% c("netVisual_diffInteraction", "netVisual_heatmap"))
        stop("unexpected return class: ", paste(class(value), collapse = ", "))
    }, finally = grDevices::dev.off())
    if (!file.exists(out) || file.info(out)$size <= 1000) stop("missing or empty PDF")
    ""
  }, error = function(e) conditionMessage(e))
  data.frame(figure = name, status = if (nzchar(detail)) "FAIL" else "PDF-WRITTEN",
             detail = detail, file = out, stringsAsFactors = FALSE)
})
summary <- do.call(rbind, results)
write.csv(summary, file.path(out.dir, "summary.csv"), row.names = FALSE)
print(summary)
if (any(summary$status != "PDF-WRITTEN")) quit(status = 1L)
