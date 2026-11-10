# Experimental local copy of metR 0.19.0's integrator. No package namespace is modified.
setwd(local({a <- commandArgs(FALSE); f <- sub('^--file=', '', grep('^--file=', a, value=TRUE)); if(length(f)) dirname(dirname(normalizePath(f))) else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {
  if (requireNamespace('renv', quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')
}
suppressPackageStartupMessages({library(methods); library(Matrix); library(cli); library(dplyr); library(tidyr); library(purrr); library(sf); library(ggplot2); library(colorspace); library(spatstat.sparse)})
source('R/SpatialCellChat_class.R'); source('R/database.R'); source('R/utilities.R')
source('R/spatial.R'); source('R/modeling.R'); source('R/analysis.R'); source('R/visualization.R')
ns <- asNamespace('metR')
original.memoised <- get('streamline.f', ns)
original.integrator <- get('_f', environment(original.memoised))
patched.integrator <- original.integrator
n.replacements <- 0L
replace.expression <- function(x) {
  if (identical(x, quote(dx + dy != 0))) {
    n.replacements <<- n.replacements + 1L
    return(quote(abs(dx) + abs(dy) != 0))
  }
  if (is.call(x)) return(as.call(lapply(as.list(x), replace.expression)))
  x
}
body(patched.integrator) <- replace.expression(body(original.integrator))
if (n.replacements != 2L) stop('Installed metR integrator differs from the 0.19.0 diagnostic; inspect it before applying this experiment.')
# The copied compute_group looks up the copied function in its own environment;
# all other symbols continue to resolve to the unmodified metR namespace.
compute.group <- get('f', environment(metR::StatStreamline$compute_group))
local.env <- new.env(parent=environment(compute.group))
local.env$streamline.f <- memoise::memoise(patched.integrator)
environment(compute.group) <- local.env
patched.stat <- ggplot2::ggproto(NULL, metR::StatStreamline, compute_group=compute.group)
coords <- as.matrix(expand.grid(x=0:4,y=0:4)); cells <- paste0('c',seq_len(nrow(coords))); rownames(coords) <- cells
expr <- Matrix::Matrix(matrix(1,nrow=1,ncol=length(cells),dimnames=list('G1',cells)),sparse=TRUE)
meta <- data.frame(label=factor(rep('A',length(cells))),row.names=cells)
chat <- createSpatialCellChat(expr,meta=meta,group.by='label',datatype='spatial',coordinates=coords,spatial.factors=list(ratio=1,tol=0.5))
out.dir <- 'tests_dev/figures/P3/streamline-patch-experiment-2026-11-11'
dir.create(out.dir,recursive=TRUE,showWarnings=FALSE)
results <- list(); warnings <- list()
capture.warnings <- function(case, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- data.frame(case=case, warning=conditionMessage(w))
    invokeRestart('muffleWarning')
  })
}
for (case in c('exact', 'near', 'opposite_diagonal')) {
  value <- switch(case, exact=1, near=.999, opposite_diagonal=-1)
  components <- matrix(rep(c(value,1),each=length(cells)),ncol=2,dimnames=list(cells,c('x','y')))
  chat@netP$cell$field$outgoing <- SparseChatArray(list(signal=as(components,'dgCMatrix')),dimnames=list(cells,c('x','y'),'signal'))
  baseline <- capture.warnings(paste0(case,'-construction'), netVisual_CommunFlow(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0))
  # Compare with the original dependency, independently of the production guard.
  baseline$layers[[2L]] <- ggplot2::ggproto(NULL, baseline$layers[[2L]])
  baseline$layers[[2L]]$stat <- metR::StatStreamline
  fixed <- baseline
  fixed$layers[[2L]] <- ggplot2::ggproto(NULL, baseline$layers[[2L]])
  fixed$layers[[2L]]$stat <- patched.stat
  baseline.lines <- capture.warnings(paste0(case,'-baseline'),ggplot2::ggplot_build(baseline)$data[[2L]])
  fixed.lines <- capture.warnings(paste0(case,'-patched'),ggplot2::ggplot_build(fixed)$data[[2L]])
  same <- isTRUE(all.equal(baseline.lines,fixed.lines,tolerance=0,check.attributes=TRUE))
  if (case != 'exact' && !same) stop('Unaffected control changed: ',case)
  if (case == 'exact' && (nrow(baseline.lines) != 0L || nrow(fixed.lines) == 0L)) stop('Expected cancellation/recovery was not observed')
  file <- file.path(out.dir,paste0(case,'-local-patch.pdf'))
  capture.warnings(paste0(case,'-pdf'),ggplot2::ggsave(file,fixed,width=6,height=5))
  results[[length(results)+1L]] <- data.frame(case=case, baseline_rows=nrow(baseline.lines),patched_rows=nrow(fixed.lines),
    unchanged_control=same, replacements=n.replacements, metR_version=as.character(packageVersion('metR')),file=file)
}
stopifnot(identical(get('streamline.f',ns),original.memoised))
summary <- do.call(rbind,results)
write.csv(summary,file.path(out.dir,'summary.csv'),row.names=FALSE)
warning.table <- if(length(warnings)) unique(do.call(rbind,warnings)) else data.frame(case=character(),warning=character())
write.csv(warning.table,file.path(out.dir,'warnings.csv'),row.names=FALSE)
print(summary)
