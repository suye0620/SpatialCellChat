# Compare the production adapter to the original metR stat on producer-derived fields.
setwd(local({a<-commandArgs(FALSE); f<-sub('^--file=','',grep('^--file=',a,value=TRUE));if(length(f))dirname(dirname(normalizePath(f)))else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {
  if (requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')
}
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(sf);library(ggplot2);library(colorspace);library(spatstat.sparse)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/spatial.R');source('R/modeling.R');source('R/analysis.R');source('R/visualization.R')
source('tests_dev/fixtures/make_simulated_chat.R');source('tests_dev/fixtures/make_visium_pp1.R');source('tests_dev/fixtures/make_slideseq_window.R')
args<-commandArgs(TRUE); if(length(args)!=1L || !args %in% c('simulated','visium','slideseq')) stop('Supply simulated, visium or slideseq')
dataset<-args[[1L]]
chat<-switch(dataset,simulated=make_simulated_chat(),visium=make_visium_chat(),slideseq=make_slideseq_chat())
signal<-dimnames(chat@netP$cell$prob)[[3L]][1L]
chat<-computeCommunField(chat,signaling.name=signal)
field.before<-chat@netP$cell$field$outgoing[[signal]]
original.memoised<-get('streamline.f',asNamespace('metR'))
out.dir<-'tests_dev/figures/P3/streamline-production-comparison-2026-11-11'
dir.create(out.dir,recursive=TRUE,showWarnings=FALSE)
warnings<-list()
capture<-function(case,expr) withCallingHandlers(expr,warning=function(w){warnings[[length(warnings)+1L]]<<-data.frame(case=case,warning=conditionMessage(w));invokeRestart('muffleWarning')})
results<-list()
for(res in c(.5,8)) {
  tag<-if(res==.5)'default' else 'res8'
  set.seed(108)
  adapted<-capture(paste(dataset,tag,'construction'),netVisual_CommunFlow(chat,signaling=signal,pattern='outgoing',res=res,point.size=.22,image.alpha=.4))
  original<-adapted; original$layers[[2L]]<-ggplot2::ggproto(NULL,adapted$layers[[2L]]);original$layers[[2L]]$stat<-metR::StatStreamline
  set.seed(108); a<-capture(paste(dataset,tag,'original-build'),ggplot2::ggplot_build(original)$data[[2L]])
  set.seed(108); b<-capture(paste(dataset,tag,'adapted-build'),ggplot2::ggplot_build(adapted)$data[[2L]])
  input<-adapted$layers[[2L]]$data
  for(kind in c('original','adapted')) {
    file<-file.path(out.dir,paste(dataset,tag,kind,sep='-')); file<-paste0(file,'.pdf')
    set.seed(108); capture(paste(dataset,tag,kind,'pdf'),ggplot2::ggsave(file,if(kind=='original')original else adapted,width=7,height=5.5))
    if(!file.exists(file)||file.info(file)$size<1000)stop('Missing PDF: ',file)
  }
  results[[length(results)+1L]]<-data.frame(dataset=dataset,signal=signal,cells=nrow(chat@meta),res=res,
    interpolated_points=nrow(input),nonzero_points=sum(input$u^2+input$v^2>0,na.rm=TRUE),
    cancelling_grid_points=sum(input$u==input$v & input$u^2+input$v^2>0,na.rm=TRUE),
    original_rows=nrow(a),adapted_rows=nrow(b),identical_drawn_layer=isTRUE(all.equal(a,b,tolerance=0,check.attributes=TRUE)))
}
stopifnot(identical(field.before,chat@netP$cell$field$outgoing[[signal]]),identical(original.memoised,get('streamline.f',asNamespace('metR'))))
write.csv(do.call(rbind,results),file.path(out.dir,paste0(dataset,'-summary.csv')),row.names=FALSE)
warning.table<-if(length(warnings))unique(do.call(rbind,warnings))else data.frame(case=character(),warning=character())
write.csv(warning.table,file.path(out.dir,paste0(dataset,'-warnings.csv')),row.names=FALSE)
print(do.call(rbind,results));print(warning.table)
