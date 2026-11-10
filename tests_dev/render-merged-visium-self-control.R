# Visium self-comparison through the actual merge producer. This is a technical
# zero-difference control, not an independent-sample biological comparison.
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));dirname(dirname(normalizePath(f)))}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {if(requireNamespace('renv',quietly=TRUE))renv::load(getwd()) else source('renv/activate.R')}
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(ggplot2);library(colorspace);library(igraph);library(cowplot);library(circlize);library(ComplexHeatmap);library(ggalluvial);library(spatstat.sparse);library(Rcpp)})
for(f in c('SpatialCellChat_class','database','utilities','spatial','modeling','merging','visualization','analysis'))source(paste0('R/',f,'.R'))
Rcpp::sourceCpp('src/SpatialChat_Rcpp.cpp', rebuild=FALSE,showOutput=FALSE)
out <- 'tests_dev/figures/merged/visium-self-control-2026-11-11';dir.create(out,recursive=TRUE,showWarnings=FALSE)
chat <- readRDS('tests_dev/test_data/visium_pp1_chat.rds')
chat <- netAnalysis_computeCentrality(chat,slot.name='netP',do.group=TRUE,degree.only=TRUE)
merged <- mergeSpatialCellChat(list(control_a=chat,control_b=chat),show.plot=FALSE)
stopifnot(isTRUE(validateSpatialCellChat(merged)),
  identical(merged@misc$.datasets$control_a$net$group$count,merged@misc$.datasets$control_b$net$group$count),
  identical(merged@misc$.datasets$control_a$net$group$weight,merged@misc$.datasets$control_b$net$group$weight),
  identical(merged@misc$.datasets$control_a$netP$centr,merged@misc$.datasets$control_b$netP$centr))
signal <- chat@netP$pathways[1L]
cell.group <- levels(chat@idents)[1L]
write.csv(data.frame(dataset='visium_PP1',copies=2,independent_samples=FALSE,cells_per_copy=nrow(chat@meta),
  groups=nlevels(chat@idents),pathways=length(chat@netP$pathways),chosen_pathway=signal),file.path(out,'fixture.csv'),row.names=FALSE)
set.seed(42)
embedding <- tryCatch({
  x <- computeNetSimilarityPairwise(merged,slot.name='netP',comparison=c(1,2),k=2)
  x <- netEmbedding(x,slot.name='netP',comparison=c(1,2),umap.method='uwot',
                    pathway.remove=character(),n_neighbors=10,n_threads=1)
  netClustering(x,slot.name='netP',comparison=c(1,2),k=2,do.plot=FALSE)
}, error=function(e)e)
require.embedding <- function() {
  if(inherits(embedding,'error')) stop('pairwise embedding setup: ',conditionMessage(embedding))
  embedding
}
cases <- list(
 compareInteractions=function()compareInteractions(merged,measure='count'),
 rankNet=function()rankNet(merged,slot.name='netP',mode='comparison',measure='count',return.data=TRUE)$gg.obj,
 netVisual_barplot=function()netVisual_barplot(merged,comparison=c(1,2),measure='count'),
 netVisual_diffInteraction=function()netVisual_diffInteraction(merged,comparison=c(1,2),measure='count'),
 netVisual_heatmap=function()netVisual_heatmap(merged,comparison=c(1,2),measure='count'),
 netVisual_bubble=function()netVisual_bubble(merged,comparison=c(1,2),signaling=signal,return.data=TRUE)$gg.obj,
 netAnalysis_diff_signalingRole_scatter=function()netAnalysis_diff_signalingRole_scatter(merged,comparison=c(1,2),do.label=FALSE,show.legend=FALSE),
 netAnalysis_signalingChanges_scatter=function()netAnalysis_signalingChanges_scatter(merged,idents.use=cell.group,comparison=c(1,2),do.label=FALSE,show.legend=FALSE),
 netVisual_embeddingPairwise=function()netVisual_embeddingPairwise(require.embedding(),slot.name='netP',comparison=c(1,2),do.label=FALSE,show.legend=FALSE),
 netVisual_embeddingPairwiseZoomIn=function()netVisual_embeddingPairwiseZoomIn(require.embedding(),slot.name='netP',comparison=c(1,2),pathway.remove=character(),do.label=FALSE,show.legend=FALSE)
)
results <- list(); warning.rows <- list()
for(name in names(cases)) {
 file <- file.path(out,paste0(name,'.pdf'));caught<-character()
 pdf(file,width=10,height=8)
 err <- tryCatch(withCallingHandlers({value<-cases[[name]]();
   if(name=='netVisual_heatmap') {m<-value@ht_list[[1L]]@matrix; stopifnot(all(m==0),identical(rownames(m),levels(chat@idents)),identical(colnames(m),levels(chat@idents)))}
   if(name=='netAnalysis_signalingChanges_scatter') stopifnot(nrow(value$data)==0)
   if(inherits(value,'ggplot'))print(value)
   else if(!name %in% c('netVisual_diffInteraction','netVisual_heatmap'))stop('unexpected return class: ',paste(class(value),collapse=', '))
   NA_character_},warning=function(w){caught<<-c(caught,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)conditionMessage(e))
 dev.off()
 if(name=='netVisual_diffInteraction' && any(grepl('zero-extent|no non-missing arguments',caught))) stop('Unexpected empty-edge warning')
 if(length(caught)) warning.rows[[length(warning.rows)+1L]]<-data.frame(figure=name,warning=caught)
 results[[length(results)+1L]]<-data.frame(figure=name,status=if(is.na(err))'DRAWN-REVIEW' else 'FAIL',
   detail=if(is.na(err))'' else err,file=file)
 print(results[[length(results)]])
}
write.csv(do.call(rbind,results),file.path(out,'summary.csv'),row.names=FALSE)
write.csv(if(length(warning.rows))do.call(rbind,warning.rows) else data.frame(figure=character(),warning=character()),
  file.path(out,'warnings.csv'),row.names=FALSE)
if(any(vapply(results,function(x)x$status=='FAIL',logical(1))))quit(status=1L)
