# Compare the public display-scale option on producer-derived fields.
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));if(length(f))dirname(dirname(normalizePath(f)))else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if(!nzchar(Sys.getenv('RENV_PROJECT'))) {if(requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')}
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(sf);library(ggplot2);library(colorspace);library(spatstat.sparse)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/spatial.R');source('R/modeling.R');source('R/analysis.R');source('R/visualization.R')
source('tests_dev/fixtures/make_simulated_chat.R');source('tests_dev/fixtures/make_visium_pp1.R');source('tests_dev/fixtures/make_slideseq_window.R')
args<-commandArgs(TRUE);if(length(args)!=1L||!args %in% c('simulated','visium','slideseq'))stop('Supply simulated, visium or slideseq')
dataset<-args[[1L]];chat<-switch(dataset,simulated=make_simulated_chat(),visium=make_visium_chat(),slideseq=make_slideseq_chat())
signal<-dimnames(chat@netP$cell$prob)[[3L]][1L];chat<-computeCommunField(chat,signaling.name=signal);field<-chat@netP$cell$field$outgoing[[signal]]
out.dir<-'tests_dev/figures/P3/arrow-display-scale-2026-11-11';dir.create(out.dir,recursive=TRUE,showWarnings=FALSE)
warnings<-list();capture<-function(case,expr)withCallingHandlers(expr,warning=function(w){warnings[[length(warnings)+1L]]<<-data.frame(case=case,warning=conditionMessage(w));invokeRestart('muffleWarning')})
plots<-list();drawn<-list()
for(scale in c(1,4)) {
  key<-as.character(scale)
  plots[[key]]<-capture(paste(dataset,key,'construction'),netVisual_CommunFieldGrid(chat,signaling=signal,pattern='outgoing',arrow.skip=0,point.size=.22,image.alpha=.4,arrow.scale=scale))
  drawn[[key]]<-capture(paste(dataset,key,'build'),ggplot_build(plots[[key]]))
  file<-file.path(out.dir,paste0(dataset,'-scale-',key,'.pdf'))
  capture(paste(dataset,key,'pdf'),ggsave(file,plots[[key]],width=7,height=5.5))
  if(!file.exists(file)||file.info(file)$size<1000)stop('Missing PDF: ',file)
}
a<-drawn[['1']]$data[[2L]];b<-drawn[['4']]$data[[2L]]
shared<-setdiff(intersect(names(a),names(b)),'mag')
same<-isTRUE(all.equal(a[,shared],b[,shared],tolerance=0,check.attributes=TRUE))
scaled<-isTRUE(all.equal(b$mag,a$mag*4,tolerance=0,check.attributes=TRUE))
stopifnot(same,scaled,identical(drawn[['1']]$data[[1L]],drawn[['4']]$data[[1L]]),
  identical(plots[['1']]$layers[[2L]]$data,plots[['4']]$layers[[2L]]$data),
  identical(field,chat@netP$cell$field$outgoing[[signal]]),
  identical(plots[['4']]$labels$caption,'Arrow length display scale: 4x'))
result<-data.frame(dataset=dataset,signal=signal,cells=nrow(chat@meta),rows=nrow(a),nonzero=sum(a$mag>0),
  median_original_cm=median(a$mag[a$mag>0]),maximum_original_cm=max(a$mag),
  median_scaled_cm=median(b$mag[b$mag>0]),maximum_scaled_cm=max(b$mag),
  other_columns_identical=same,magnitude_times_four=scaled)
write.csv(result,file.path(out.dir,paste0(dataset,'-summary.csv')),row.names=FALSE)
write.csv(if(length(warnings))unique(do.call(rbind,warnings))else data.frame(case=character(),warning=character()),file.path(out.dir,paste0(dataset,'-warnings.csv')),row.names=FALSE)
print(result)
