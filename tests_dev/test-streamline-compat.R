# Regression tests for the production metR streamline compatibility adapter.
setwd(local({a <- commandArgs(FALSE); f <- sub('^--file=', '', grep('^--file=', a, value=TRUE)); if(length(f)) dirname(dirname(normalizePath(f))) else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {
  if (requireNamespace('renv', quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')
}
suppressPackageStartupMessages({library(methods); library(Matrix); library(cli); library(dplyr); library(tidyr); library(purrr); library(sf); library(ggplot2); library(colorspace); library(spatstat.sparse)})
source('R/SpatialCellChat_class.R'); source('R/database.R'); source('R/utilities.R')
source('R/spatial.R'); source('R/modeling.R'); source('R/analysis.R'); source('R/visualization.R')
expect_error <- function(expr, text) {
  result <- tryCatch({force(expr); NULL}, error=function(e) conditionMessage(e))
  if (is.null(result) || !grepl(text, result, fixed=TRUE)) stop('Expected error containing: ', text, '; got: ', result)
}
ns <- asNamespace('metR')
old.integrator <- get('streamline.f', ns, inherits=FALSE)
old.stat <- metR::StatStreamline
old.body <- body(get('_f', environment(old.integrator), inherits=FALSE))
# A probe must not create or advance the user's random state.
set.seed(421); before <- .Random.seed
invisible(.sc_streamline_probe(old.stat, -.999))
stopifnot(identical(before, .Random.seed))
rm('.Random.seed', envir=.GlobalEnv)
invisible(.sc_streamline_probe(old.stat, -.999))
stopifnot(!exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE))
# Verify each version-selection and refusal branch with controlled probes.
nonempty <- data.frame(x=1)
empty <- nonempty[FALSE,,drop=FALSE]
fixed.probe <- function(stat, dy) nonempty
stopifnot(identical(.sc_resolve_streamline_stat(old.stat, package_version('99.0'), ns, fixed.probe), old.stat))
affected.probe <- function(stat, dy) if (dy == -1) empty else nonempty
expect_error(.sc_resolve_streamline_stat(old.stat, package_version('99.0'), ns, affected.probe), 'has not been verified')
expect_error(.sc_resolve_streamline_stat(old.stat, package_version('0.19.0'), new.env(), affected.probe), 'cannot be safely corrected')
mismatch.ns <- new.env(); mismatch.env <- new.env(); mismatch.env$`_f` <- function(dx,dy) dx+dy
mismatch.memoised <- function(...) NULL; environment(mismatch.memoised) <- mismatch.env
mismatch.ns$streamline.f <- mismatch.memoised
expect_error(.sc_resolve_streamline_stat(old.stat, package_version('0.19.0'), mismatch.ns, affected.probe), 'expected exactly two zero-vector predicates')
expect_error(.sc_resolve_streamline_stat(old.stat, package_version('0.19.0'), ns, affected.probe), 'did not pass its compatibility checks')
expect_error(.sc_resolve_streamline_stat(old.stat, package_version('0.19.0'), ns, function(stat, dy) empty), 'control field')
custom.stat <- ggplot2::ggproto(NULL, old.stat)
custom.layer <- metR::geom_streamline(stat=custom.stat)
stopifnot(identical(.sc_guard_streamline_layer(custom.layer)$stat, custom.stat))
stat <- .sc_streamline_compat_stat()
stopifnot(identical(stat, .sc_streamline_compat_stat()))
# A genuinely corrected stat must pass the real behavioral probe unchanged.
stopifnot(identical(.sc_resolve_streamline_stat(stat, package_version('99.0'), new.env()), stat))
# Compare real metR and production adapter, with reproducible seeds for jitter.
make.field <- function(case) {
  d <- expand.grid(x=-3:3,y=-3:3)
  if(case=='zero') {d$dx<-0;d$dy<-0}
  if(case=='axis') {d$dx<-1;d$dy<-0}
  if(case=='cancel') {d$dx<-1;d$dy<--1}
  if(case=='near') {d$dx<-1;d$dy<--.999}
  if(case=='curved') {d$dx<--.3*d$y;d$dy<-.3*d$x+.17}
  if(case %in% c('heterogeneous','boundary','periodic','jitter')) {d$dx<-1+.04*d$x;d$dy<-.4+.03*d$y}
  d
}
run.build <- function(p) {
  warnings <- character(); set.seed(103)
  value <- tryCatch(withCallingHandlers(ggplot_build(p)$data[[1L]],warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  list(value=value,warnings=warnings)
}
cases <- c('zero','axis','cancel','near','curved','heterogeneous','boundary','periodic','jitter')
results <- lapply(cases, function(case) {
  d <- make.field(case)
  start <- if(case=='boundary') data.frame(x=2.99,y=.2) else data.frame(x=.7,y=.4)
  params <- list(mapping=aes(x=x,y=y,dx=dx,dy=dy),L=2,res=8,dt=.05,S=12,min.L=0,start=start,jitter=if(case=='jitter') 0.3 else 0)
  if(case=='periodic') params$xwrap<-c(-3,3)
  p <- ggplot(d)+do.call(metR::geom_streamline,params)
  fixed <- p; fixed$layers[[1L]] <- ggplot2::ggproto(NULL,p$layers[[1L]]); fixed$layers[[1L]]$stat <- stat
  a <- run.build(p); b <- run.build(fixed)
  same <- isTRUE(all.equal(a$value,b$value,tolerance=0,check.attributes=TRUE)) && identical(a$warnings,b$warnings)
  if(case!='cancel' && !same) stop('Control changed: ',case)
  if(case=='cancel' && (inherits(b$value,'error') || !nrow(b$value))) stop('Cancellation was not recovered')
  data.frame(case=case,baseline_rows=if(inherits(a$value,'error')) NA_integer_ else nrow(a$value),
    adapted_rows=if(inherits(b$value,'error')) NA_integer_ else nrow(b$value),unchanged=same,
    baseline_warning=paste(a$warnings,collapse=' | '),adapted_warning=paste(b$warnings,collapse=' | '))
})
# Distinguish exactly zero from small nonzero or invalid values.
stopifnot(.sc_is_zero_streamline_data(data.frame(u=c(0,0),v=c(0,0))),
  !.sc_is_zero_streamline_data(data.frame(u=c(0,1e-300),v=c(0,0))),
  !.sc_is_zero_streamline_data(data.frame(u=c(0,NA_real_),v=c(0,0))),
  !.sc_is_zero_streamline_data(data.frame(u=c(0,Inf),v=c(0,0))),
  !.sc_is_zero_streamline_data(data.frame(u=numeric(),v=numeric())),
  !.sc_is_zero_streamline_data(data.frame(u=0)))
# Exercise the actual public entry point and a serialized/reloaded plot.
coords <- as.matrix(expand.grid(x=0:4,y=0:4)); cells <- paste0('c',seq_len(nrow(coords))); rownames(coords) <- cells
expr <- Matrix::Matrix(matrix(1,nrow=1,ncol=length(cells),dimnames=list('G1',cells)),sparse=TRUE)
meta <- data.frame(label=factor(rep('A',length(cells))),row.names=cells)
chat <- createSpatialCellChat(expr,meta=meta,group.by='label',datatype='spatial',coordinates=coords,spatial.factors=list(ratio=1,tol=.5))
components <- matrix(rep(c(1,1),each=length(cells)),ncol=2,dimnames=list(cells,c('x','y')))
chat@netP$cell$field$outgoing <- SparseChatArray(list(signal=as(components,'dgCMatrix')),dimnames=list(cells,c('x','y'),'signal'))
p <- netVisual_CommunFlow(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0)
custom.plot <- netVisual_CommunFlow(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0,stat=custom.stat)
stopifnot(identical(custom.plot$layers[[2L]]$stat,custom.stat))
zero.chat <- chat
zero.components <- Matrix::Matrix(components * 0,sparse=TRUE)
zero.chat@netP$cell$field$outgoing <- SparseChatArray(list(signal=as(zero.components,'dgCMatrix')),dimnames=list(cells,c('x','y'),'signal'))
zero.before <- zero.chat
zero.plot <- suppressWarnings(netVisual_CommunFlow(zero.chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0))
zero.warnings <- character()
zero.built <- withCallingHandlers(ggplot2::ggplot_build(zero.plot),warning=function(w){zero.warnings<<-c(zero.warnings,conditionMessage(w));invokeRestart('muffleWarning')})
stopifnot(identical(zero.chat,zero.before),nrow(zero.built$data[[1L]])==length(cells),
  nrow(zero.built$data[[2L]])==0L,nrow(zero.plot$layers[[2L]]$data)==0L,
  !any(grepl('stat_streamline|integrat',zero.warnings,ignore.case=TRUE)))
zero.pdf <- tempfile(fileext='.pdf'); ggplot2::ggsave(zero.pdf,zero.plot,width=6,height=5)
stopifnot(file.info(zero.pdf)$size>0); unlink(zero.pdf)
zero.custom <- suppressWarnings(netVisual_CommunFlow(zero.chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0,stat=custom.stat))
stopifnot(identical(zero.custom$layers[[2L]]$stat,custom.stat),nrow(zero.custom$layers[[2L]]$data)>0L,
  all(zero.custom$layers[[2L]]$data$u==0),all(zero.custom$layers[[2L]]$data$v==0))
built <- ggplot2::ggplot_build(p)$data[[2L]]
stopifnot(nrow(built)>0)
file <- tempfile(fileext='.rds'); pdf.file <- tempfile(fileext='.pdf')
saveRDS(p,file); restored <- readRDS(file)
rebuilt <- ggplot2::ggplot_build(restored)$data[[2L]]
stopifnot(isTRUE(all.equal(built,rebuilt,tolerance=0,check.attributes=TRUE)))
ggplot2::ggsave(pdf.file,restored,width=6,height=5)
stopifnot(file.info(pdf.file)$size>0)
child <- tempfile(fileext='.R')
writeLines(c('suppressPackageStartupMessages(library(ggplot2))',
  'p <- readRDS(commandArgs(TRUE)[1L])',
  'd <- ggplot2::ggplot_build(p)$data[[2L]]',
  'if (!nrow(d)) stop("No streamlines after a fresh-process reload")',
  'ggplot2::ggsave(commandArgs(TRUE)[2L], p, width=6, height=5)'), child)
child.status <- system2(file.path(R.home('bin'),'Rscript.exe'),c(shQuote(child),shQuote(file),shQuote(pdf.file)))
stopifnot(identical(child.status,0L),file.info(pdf.file)$size>0)
unlink(c(file,pdf.file,child))
stopifnot(identical(get('streamline.f',ns,inherits=FALSE),old.integrator), identical(metR::StatStreamline,old.stat), identical(body(get('_f',environment(old.integrator),inherits=FALSE)),old.body))
desc <- read.dcf('DESCRIPTION'); stopifnot(grepl('memoise',desc[1,'Imports'],fixed=TRUE))
lock <- jsonlite::fromJSON('renv.lock',simplifyVector=FALSE); stopifnot(!is.null(lock$Packages$memoise))
print(do.call(rbind,results)); cat('Production adapter branch, state, public entry and serialization checks passed.\n')
