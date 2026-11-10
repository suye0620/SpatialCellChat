args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, args[1L] %in% c('simulated','visium','slideseq'))
dataset <- args[1L]
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));dirname(dirname(normalizePath(f)))}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) { if (requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R') }
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(ggplot2);library(colorspace);library(circlize);library(ComplexHeatmap);library(grid);library(spatstat.sparse)})
for (f in c('SpatialCellChat_class','database','utilities','spatial','modeling','analysis','visualization')) source(paste0('R/',f,'.R'))
fixture <- switch(dataset, simulated='tests_dev/test_data/simulated_chat.rds', visium='tests_dev/test_data/visium_pp1_chat.rds',slideseq='tests_dev/test_data/slideseq_full_chat.rds')
chat <- readRDS(fixture);signal <- chat@netP$pathways[1L]
out <- file.path('tests_dev/figures/P2/chord-production-outside-2026-11-11',dataset)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Intercept only chordDiagram calls; capture evaluated inputs after drawing.
current <- NULL
chordDiagram <- function(df, ...) {
  ans <- circlize::chordDiagram(df, ...)
  current <<- list(df=df,args=list(...))
  ans
}
old.expr <- parse(text=system2('git',c('show','f643c28:R/visualization.R'),stdout=TRUE))
old.def <- Filter(function(e)is.call(e) && identical(e[[1]],as.name('<-')) &&
  identical(e[[2]],as.name('netVisual_chord_gene')),old.expr)
stopifnot(length(old.def)==1L)
eval(old.def[[1]]);old.fun <- netVisual_chord_gene
source('R/visualization.R')
# Instrument a copy of the actual outside helper to inspect its rendered
# paths and labels without changing their drawing operations.
probe <- new.env(parent=emptyenv())
record.outside.lines <- function(x,y,col,...) {
  probe$paths[[length(probe$paths)+1L]] <- list(x=x,y=y,col=col)
  graphics::lines(x,y,col=col,...)
}
record.outside.text <- function(x,y,labels,...) {
  probe$text[[length(probe$text)+1L]] <- list(x=x,y=y,label=labels)
  graphics::text(x,y,labels=labels,...)
}
rewrite <- function(x) {
  if(is.call(x) && is.call(x[[1]]) && identical(x[[1]][[1]],as.name('::')) &&
     identical(x[[1]][[2]],as.name('graphics')) && as.character(x[[1]][[3]]) %in% c('lines','text')) {
    x[[1]] <- as.name(if(as.character(x[[1]][[3]])=='lines') 'record.outside.lines' else 'record.outside.text')
  }
  if(is.call(x)) return(as.call(lapply(as.list(x),rewrite)))
  x
}
helper <- .sc_chord_gene_outside_labels
expr <- as.list(body(helper))
loop <- which(vapply(expr,function(x) is.call(x) && identical(x[[1]],as.name('for')) &&
  grepl('seq_along(sectors)',paste(deparse(x),collapse=' '),fixed=TRUE),logical(1)))
stopifnot(length(loop)==1L)
expr <- append(expr,list(quote(probe$meta <- list(sectors=sectors,x=x,y=y,radius=radius,label=label,label.y=label.y))),after=loop-1L)
body(helper) <- rewrite(as.call(expr))
.sc_chord_gene_outside_labels <- helper
cases <- data.frame(kind=c('old','radial','outside-large','outside-small','outside-small-legend-overlap'),
                    width=c(24,24,24,9,9),height=c(16,16,16,7,7),
                    outside=c(FALSE,FALSE,TRUE,TRUE,TRUE))
inputs <- list();results <- list();warnings <- list()
for(i in seq_len(nrow(cases))) {
  case <- cases[i,]; file <- file.path(out,paste0(case$kind,'.pdf'))
  grDevices::pdf(file,width=case$width,height=case$height)
  fun <- if(case$kind=='old') old.fun else netVisual_chord_gene
  callargs <- list(object=chat,signaling=signal,show.legend=TRUE,lab.cex=.7,small.gap=.25)
  if(case$outside) callargs$label.layout <- 'outside'
  if(case$kind=='outside-small-legend-overlap') callargs$legend.pos.x <- 30
  current <- NULL;captured <- character();probe$paths <- list();probe$text <- list();probe$meta <- NULL
  err <- tryCatch(withCallingHandlers({result<-do.call(fun,callargs);stopifnot(inherits(result,'recordedplot'));NA_character_},
     warning=function(w){captured<<-c(captured,conditionMessage(w));invokeRestart('muffleWarning')}),
     error=function(e)conditionMessage(e))
  grDevices::dev.off();stopifnot(length(circlize::get.all.sector.index())==0L)
  inputs[[case$kind]] <- current
  if (case$outside && is.na(err)) {
    meta <- probe$meta; stopifnot(length(meta$sectors)==length(probe$paths),length(probe$text)==length(meta$sectors),!anyDuplicated(meta$sectors))
    for (k in seq_along(meta$sectors)) {
      path <- probe$paths[[k]]; text <- probe$text[[k]]
      stopifnot(identical(text$label,meta$label[k]),
        isTRUE(all.equal(unname(path$x[1]),unname(meta$x[k]*meta$radius[k]))),
        isTRUE(all.equal(unname(path$y[1]),unname(meta$y[k]*meta$radius[k]))),
        isTRUE(all.equal(unname(tail(path$y,1)),unname(text$y))),
        identical(unname(path$col),unname(current$args$grid.col[meta$sectors[k]])))
    }
    paths <- do.call(rbind,lapply(seq_along(meta$sectors),function(k) data.frame(
      sector=meta$sectors[k],vertex=seq_along(probe$paths[[k]]$x),x=probe$paths[[k]]$x,y=probe$paths[[k]]$y)))
    write.csv(paths,file.path(out,paste0(case$kind,'-paths.csv')),row.names=FALSE)
  }
  results[[i]] <- data.frame(dataset=dataset,case=case$kind,signal=signal,status=if(is.na(err))'DRAWN-REVIEW' else 'ERROR',
     error=if(is.na(err))'' else err,warnings=paste(captured,collapse=' | '),file=file)
  print(results[[i]])
}
a<-inputs$old;b<-inputs$radial;stopifnot(identical(a,b))
a$args$preAllocateTracks<-NULL
for(kind in c('outside-large','outside-small','outside-small-legend-overlap')) {b<-inputs[[kind]];b$args$preAllocateTracks<-NULL;stopifnot(identical(a,b))}
stopifnot(all(vapply(results,function(x)x$status=='DRAWN-REVIEW',logical(1))))
# These conditions specifically reproduce the measured Visium failure.
if(dataset=='visium') {
 stopifnot(grepl('exceed available space',results[[4]]$warnings,fixed=TRUE),
           grepl('overlap the legend',results[[5]]$warnings,fixed=TRUE),
           !grepl('exceed available space',results[[3]]$warnings,fixed=TRUE))
}
# Invalid layout and invalid fontsize must fail without leaving circlize state.
for(params in list(list(label.layout='bad'),list(label.layout='outside',lab.cex=0))) {
  failed <- inherits(try(do.call(netVisual_chord_gene,c(list(object=chat,signaling=signal),params)),silent=TRUE),'try-error')
  stopifnot(failed,length(circlize::get.all.sector.index())==0L)
}
if(dataset=='simulated') {
  # Pixel equality checks the default rendering, not only its chord arguments.
  images <- file.path(tempdir(),c('old-chord.png','new-chord.png'))
  for (j in 1:2) {
    grDevices::png(images[j], width=1200, height=800)
    f <- if(j==1) old.fun else netVisual_chord_gene
    f(chat, 'net', NULL, signal, show.legend=FALSE, lab.cex=.7, small.gap=.25, link.sort=TRUE)
    grDevices::dev.off()
  }
  stopifnot(identical(unname(tools::md5sum(images)[1]),unname(tools::md5sum(images)[2])))
  empty <- data.frame(source=character(),target=character(),ligand=character(),receptor=character(),prob=numeric())
  stopifnot(inherits(try(netVisual_chord_gene(chat,net=empty,label.layout='outside'),silent=TRUE),'try-error'))
  long.net <- .sc_viz_group_edges(chat, .05)
  long.net$ligand <- paste0(long.net$ligand, strrep('Z', 100))
  observed <- character()
  grDevices::pdf(file.path(tempdir(),'long-label.pdf'),width=9,height=7)
  withCallingHandlers(netVisual_chord_gene(chat, net=long.net, label.layout='outside', show.legend=FALSE),
    warning=function(w) {observed <<- c(observed,conditionMessage(w));invokeRestart('muffleWarning')})
  grDevices::dev.off()
  stopifnot(any(grepl('exceed available space',observed,fixed=TRUE)))
  grDevices::pdf(file.path(tempdir(),'legend-error.pdf'))
  failed <- inherits(try(netVisual_chord_gene(chat,signaling=signal,label.layout='outside',legend.pos.x=NA_real_),silent=TRUE),'try-error')
  stopifnot(failed,length(circlize::get.all.sector.index())==0L)
  grDevices::dev.off()
}
write.csv(do.call(rbind,results),file.path(out,'summary.csv'),row.names=FALSE)
saveRDS(inputs,file.path(out,'inputs.rds'))
cat('PASS: chord input equivalence, recorded plots, warning and state checks\n')
