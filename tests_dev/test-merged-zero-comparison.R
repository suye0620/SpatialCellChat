# Verify zero-difference rendering without changing legacy directional filtering.
source('tests_dev/test-merge-spatialcellchat.R',local=FALSE)
zero <- merged
for(slot in c('net','netP')) {
  zero@misc$.datasets$sample2[[slot]]$group <- zero@misc$.datasets$sample1[[slot]]$group
  zero@misc$.datasets$sample2[[slot]]$centr <- zero@misc$.datasets$sample1[[slot]]$centr
}
for(measure in c('count','weight')) {
  grDevices::pdf(file.path(tempdir(),paste0('zero-',measure,'.pdf')))
  heat <- netVisual_heatmap(zero,measure=measure)
  grDevices::dev.off()
  mat <- heat@ht_list[[1L]]@matrix
  stopifnot(all(mat==0),identical(rownames(mat),levels(zero@idents)),identical(colnames(mat),levels(zero@idents)))
}
# The old single-object error and an all-isolates comparison remain errors.
pdf(file.path(tempdir(),'isolated.pdf'))
empty.error <- inherits(try(netVisual_heatmap(zero,measure='count',remove.isolate=TRUE),silent=TRUE),'try-error')
dev.off();stopifnot(empty.error)
single.zero <- one;single.zero@net$group$count[,] <- 0
pdf(file.path(tempdir(),'single-zero.pdf'))
single.error <- inherits(try(netVisual_heatmap(single.zero,measure='count'),silent=TRUE),'try-error')
dev.off();stopifnot(single.error)

has_message <- function(p) {
  layers <- ggplot2::ggplot_build(p)$data
  any(vapply(layers, function(layer) 'label' %in% names(layer) &&
    any(layer$label == 'No pathways passed the current change filter'), logical(1)))
}
for(label in c(FALSE,TRUE)) {
 p <- netAnalysis_signalingChanges_scatter(zero,idents.use='A',do.label=label)
 stopifnot(inherits(p,'ggplot'),nrow(p$data)==0,
           has_message(p))
 ggplot2::ggplot_build(p)
}
# Opposite nonzero directions can be removed by the legacy row-sum filter.
# The empty message must describe the filter instead of claiming no change.
cancel <- zero
for(dataset in c('sample1','sample2')) {
  cancel@misc$.datasets[[dataset]]$netP$centr['outdeg','A','P1'] <- 3
  cancel@misc$.datasets[[dataset]]$netP$centr['indeg','A','P1'] <- 3
}
cancel@misc$.datasets$sample2$netP$centr['outdeg','A','P1'] <- 5
cancel@misc$.datasets$sample2$netP$centr['indeg','A','P1'] <- 1
p <- netAnalysis_signalingChanges_scatter(cancel,idents.use='A',signaling='P1',do.label=TRUE)
stopifnot(nrow(p$data)==0, has_message(p))

# No edges should retain the nodes and avoid warnings caused by empty edge data.
warnings <- character();pdf(file.path(tempdir(),'network-zero.pdf'))
network <- withCallingHandlers(netVisual_diffInteraction(zero,measure='count'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
dev.off();stopifnot(inherits(network,'recordedplot'),
  !any(grepl('zero-extent|no non-missing arguments',warnings)))
cat('PASS: zero comparison heatmaps, filtered scatter, cancellation and edgeless network\n')
