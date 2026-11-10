# Boundary study of the isolated chord-label prototype. Does not alter production.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, args[1L] %in% c('simulated','visium','slideseq'))
dataset <- args[1L]
setwd(local({a<-commandArgs(FALSE);f<-sub('^--file=','',grep('^--file=',a,value=TRUE));dirname(dirname(normalizePath(f)))}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) { if (requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R') }
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(ggplot2);library(colorspace);library(circlize);library(ComplexHeatmap);library(grid);library(spatstat.sparse)})
for (f in c('SpatialCellChat_class','database','utilities','spatial','modeling','analysis','visualization')) source(paste0('R/',f,'.R'))
# Import just the experiment's function definitions, not its executable main body.
prototype <- parse('tests_dev/diagnose-chord-external-unfold.R')
for (expr in prototype) if (is.call(expr) && identical(expr[[1]],as.name('<-')) &&
    is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% c('plot.labels','capture.chord','make.variant')) eval(expr)
fixture <- switch(dataset, simulated='tests_dev/test_data/simulated_chat.rds',
  visium='tests_dev/test_data/visium_pp1_chat.rds',slideseq='tests_dev/test_data/slideseq_full_chat.rds')
chat <- readRDS(fixture)
signal <- chat@netP$pathways[1L]
base.out <- file.path('tests_dev/figures/P2/chord-external-boundaries-2026-11-11',dataset)
dir.create(base.out,recursive=TRUE,showWarnings=FALSE)
cases <- data.frame(name=c('large-no-legend','small-no-legend','large-legend','small-legend'),
                    width=c(24,9,24,9),height=c(16,7,16,7),legend=c(FALSE,FALSE,TRUE,TRUE))
results <- list();warnings <- list()
for (i in seq_len(nrow(cases))) {
  case <- cases[i,]
  out <- file.path(base.out,case$name);dir.create(out,recursive=TRUE,showWarnings=FALSE)
  current.kind <- 'external'; file <- file.path(out,'external.pdf')
  grDevices::pdf(file,width=case$width,height=case$height)
  error <- tryCatch(withCallingHandlers({make.variant(TRUE)(chat,signaling=signal,show.legend=case$legend,
      lab.cex=.7,small.gap=.25); NA_character_},warning=function(w) {
        warnings[[length(warnings)+1L]] <<- data.frame(case=case$name,warning=conditionMessage(w))
        invokeRestart('muffleWarning')
      }),error=function(e) conditionMessage(e))
  grDevices::dev.off();circlize::circos.clear()
  labels <- if(file.exists(file.path(out,'labels.rds'))) readRDS(file.path(out,'labels.rds')) else NULL
  overlaps <- NA_integer_;outside <- NA_integer_; n <- NA_integer_
  if (!is.null(labels)) {
    p <- labels$positions;n <- nrow(p); overlaps <- 0L
    for (side in c('left','right')) {
      x <- p[p$side==side,,drop=FALSE];x <- x[order(x$label.y),,drop=FALSE]
      if(nrow(x)>1) overlaps <- overlaps + sum(diff(x$label.y) <= (head(x$label.height,-1)+tail(x$label.height,-1))/2)
    }
    outside <- sum(1.22+p$label.width>2.1 | abs(p$label.y)+p$label.height/2>1.2)
    write.csv(p,file.path(out,'positions.csv'),row.names=FALSE)
    write.csv(labels$mapping,file.path(out,'mapping.csv'),row.names=FALSE)
  }
  results[[i]] <- data.frame(dataset=dataset,case=case$name,signal=signal,
    status=if(is.na(error)) 'DRAWN-REVIEW' else 'ERROR',error=if(is.na(error)) '' else error,
    sectors=n,adjacent_text_overlaps=overlaps,outside_canvas=outside,file=file)
  print(results[[i]])
}
write.csv(do.call(rbind,results),file.path(base.out,'summary.csv'),row.names=FALSE)
write.csv(if(length(warnings))do.call(rbind,warnings)else data.frame(case=character(),warning=character()),
          file.path(base.out,'warnings.csv'),row.names=FALSE)
