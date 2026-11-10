# Reproduce the installed metR streamline cancellation in the public flow plot.
setwd(local({a <- commandArgs(FALSE); f <- sub('^--file=', '', grep('^--file=', a, value=TRUE)); if(length(f)) dirname(dirname(normalizePath(f))) else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY = 'renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {
  if (requireNamespace('renv', quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')
}
suppressPackageStartupMessages({library(methods); library(Matrix); library(cli); library(dplyr); library(tidyr); library(purrr); library(sf); library(ggplot2); library(colorspace); library(spatstat.sparse)})
source('R/SpatialCellChat_class.R'); source('R/database.R'); source('R/utilities.R')
source('R/spatial.R'); source('R/modeling.R'); source('R/analysis.R'); source('R/visualization.R')
coords <- as.matrix(expand.grid(x=0:4, y=0:4)); cells <- paste0('c', seq_len(nrow(coords))); rownames(coords) <- cells
expr <- Matrix::Matrix(matrix(1, nrow=1, ncol=length(cells), dimnames=list('G1', cells)), sparse=TRUE)
meta <- data.frame(label=factor(rep('A', length(cells))), row.names=cells)
chat <- createSpatialCellChat(expr, meta=meta, group.by='label', datatype='spatial', coordinates=coords, spatial.factors=list(ratio=1, tol=0.5))
out.dir <- 'tests_dev/figures/P3/streamline-cancellation-2026-11-11'
dir.create(out.dir, recursive=TRUE, showWarnings=FALSE)
results <- list(); warnings <- list()
capture.warnings <- function(case, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- data.frame(case=case, warning=conditionMessage(w))
    invokeRestart('muffleWarning')
  })
}
for (value in c(1, .999)) {
  components <- matrix(rep(c(value, 1), each=length(cells)), ncol=2, dimnames=list(cells, c('x','y')))
  chat@netP$cell$field$outgoing <- SparseChatArray(list(signal=as(components,'dgCMatrix')), dimnames=list(cells, c('x','y'), 'signal'))
  case <- if (value == 1) 'exact-cancellation' else 'near-cancellation'
  plot <- capture.warnings(paste0(case, '-construction'), netVisual_CommunFlow(chat, signaling='signal', pattern='outgoing', cellsize=1,
    grid.resolution=1, L=2, min.L=0, res=8, start=data.frame(x=2,y=-2), jitter=0))
  # Keep the historical baseline explicit after the production adapter is enabled.
  plot$layers[[2L]] <- ggplot2::ggproto(NULL, plot$layers[[2L]])
  plot$layers[[2L]]$stat <- metR::StatStreamline
  uv <- plot$layers[[2]]$data
  lines <- capture.warnings(paste0(case, '-build'), ggplot2::ggplot_build(plot)$data[[2]])
  file <- file.path(out.dir, paste0(case, '.pdf'))
  capture.warnings(paste0(case, '-pdf'), ggplot2::ggsave(file, plot, width=6, height=5))
  results[[length(results)+1L]] <- data.frame(case=case, metR_version=as.character(packageVersion('metR')),
    input_x=value, input_y=1, uv_nonzero=sum(uv$u^2+uv$v^2>0, na.rm=TRUE),
    uv_exact_cancellation=sum(uv$u-uv$v==0 & uv$u^2+uv$v^2>0, na.rm=TRUE),
    line_rows=nrow(lines), distinct_positions=if (nrow(lines)) nrow(unique(lines[,c('x','y')])) else 0,
    file=file)
}
summary <- do.call(rbind, results)
write.csv(summary, file.path(out.dir, 'summary.csv'), row.names=FALSE)
warning.table <- if (length(warnings)) unique(do.call(rbind, warnings)) else data.frame(case=character(), warning=character())
write.csv(warning.table, file.path(out.dir, 'warnings.csv'), row.names=FALSE)
print(summary)
if (summary$line_rows[2] == 0) stop('near-cancellation control did not produce a streamline')
