# Experiment only: test circlize's label spreading on the full Visium pathway.
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));if(length(f))dirname(dirname(normalizePath(f)))else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if(!nzchar(Sys.getenv('RENV_PROJECT'))) {if(requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')}
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(ggplot2);library(colorspace);library(circlize);library(ComplexHeatmap);library(grid);library(spatstat.sparse)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/spatial.R');source('R/modeling.R');source('R/analysis.R');source('R/visualization.R');source('tests_dev/fixtures/make_visium_pp1.R')
chat<-make_visium_chat();signal<-chat@netP$pathways[1L]
spread.labels <- function(lab.cex) {
  sectors<-circlize::get.all.sector.index()
  positions<-vapply(sectors,function(s)mean(circlize::get.cell.meta.data('xlim',sector.index=s,track.index=1L)),numeric(1))
  circlize::circos.labels(sectors,positions,labels=trimws(sectors),cex=lab.cex,
    side='outside',connection_height=circlize::mm_h(2),padding=.2)
  data.frame(sector=sectors,x=positions)
}
experiment<-netVisual_chord_gene
count<-0L
rewrite<-function(x) {
 if(is.call(x)&&identical(x[[1L]],as.name('circos.track'))) {count<<-count+1L;return(quote(spread.labels(lab.cex)))}
 if(is.call(x))return(as.call(lapply(as.list(x),rewrite)))
 x
}
body(experiment)<-rewrite(body(experiment));stopifnot(count==1L)
out.dir<-'tests_dev/figures/P2/visium/chord-label-experiment-2026-11-11';dir.create(out.dir,recursive=TRUE,showWarnings=FALSE)
results<-list();warns<-list()
for(kind in c('original','spread')) {
 file<-file.path(out.dir,paste0(kind,'.pdf'));grDevices::pdf(file,width=20,height=20)
 outcome<-tryCatch(withCallingHandlers({if(kind=='original') netVisual_chord_gene(chat,signaling=signal,show.legend=FALSE,lab.cex=.4,small.gap=.25) else experiment(chat,signaling=signal,show.legend=FALSE,lab.cex=.4,small.gap=.25);'ok'},warning=function(w){warns[[length(warns)+1L]]<<-data.frame(kind=kind,warning=conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)conditionMessage(e))
 grDevices::dev.off();circlize::circos.clear()
 results[[length(results)+1L]]<-data.frame(kind=kind,signal=signal,result=outcome,file=file)
}
write.csv(do.call(rbind,results),file.path(out.dir,'summary.csv'),row.names=FALSE)
write.csv(if(length(warns))do.call(rbind,warns)else data.frame(kind=character(),warning=character()),file.path(out.dir,'warnings.csv'),row.names=FALSE)
print(do.call(rbind,results));print(warns)
