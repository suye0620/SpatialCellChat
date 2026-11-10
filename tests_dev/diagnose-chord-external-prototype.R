# Isolated prototype: external labels using sector mid-angles and ordered leader lines.
setwd(local({a <- commandArgs(FALSE); f <- sub('^--file=', '', grep('^--file=', a, value=TRUE)); if(length(f)) dirname(dirname(normalizePath(f))) else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) { if (requireNamespace('renv', quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R') }
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(ggplot2);library(colorspace);library(circlize);library(ComplexHeatmap);library(grid);library(spatstat.sparse)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/spatial.R');source('R/modeling.R');source('R/analysis.R');source('R/visualization.R');source('tests_dev/fixtures/make_visium_pp1.R')
out <- 'tests_dev/figures/P2/visium/chord-external-prototype-2026-11-11'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
chat <- make_visium_chat(); signal <- chat@netP$pathways[1L]
plot.labels <- function(lab.cex, grid.col, order.source, order.target) {
  sectors <- circlize::get.all.sector.index()
  mid <- vapply(sectors, function(s) mean(circlize::get.cell.meta.data('xplot', sector.index=s, track.index=1L)), numeric(1)) * pi/180
  radius <- vapply(sectors, function(s) circlize::get.cell.meta.data('cell.top.radius', sector.index=s, track.index=2L), numeric(1))
  xy <- data.frame(sector=sectors, angle=mid, radius=radius, x=cos(mid), y=sin(mid))
  xy$side <- ifelse(xy$x >= 0, 'right', 'left')
  xy$label.y <- NA_real_
  xy$label.height <- graphics::strheight(trimws(sectors), cex=lab.cex)
  xy$label.width <- graphics::strwidth(trimws(sectors), cex=lab.cex)
  for (side in c('left', 'right')) {
    ii <- which(xy$side == side); ii <- ii[order(xy$y[ii])]
    if (length(ii)) xy$label.y[ii] <- if(length(ii)==1L) 0 else seq(-1.08, 1.08, length.out=length(ii))
  }
  for (i in seq_len(nrow(xy))) {
    side <- if (xy$side[i]=='right') 1 else -1
    graphics::lines(c(xy$x[i]*xy$radius[i], xy$x[i]*1.01, side*1.12, side*1.2), c(xy$y[i]*xy$radius[i], xy$y[i]*1.01, xy$label.y[i], xy$label.y[i]), col=grid.col[xy$sector[i]], lwd=.45, xpd=NA)
    graphics::text(side*1.22, xy$label.y[i], labels=trimws(xy$sector[i]), adj=c(if(side>0) 0 else 1,.5), cex=lab.cex, xpd=NA)
  }
  source.map <- data.frame(sector=order.source$ligand, group=as.character(order.source$source), role='source')
  target.map <- data.frame(sector=order.target$receptor, group=as.character(order.target$target), role='target')
  map <- rbind(source.map, target.map)
  stopifnot(!anyDuplicated(sectors), !anyDuplicated(map$sector), setequal(sectors,map$sector),
            all(!is.na(grid.col[sectors])))
  saveRDS(list(positions=xy, mapping=map, colors=grid.col), file.path(out, 'labels.rds'))
}
capture.chord <- function(df, ...) {
  result <- circlize::chordDiagram(df, ...)
  args <- list(...); saveRDS(list(df=df, args=args), file.path(out, paste0(current.kind, '-chord-input.rds')))
  invisible(result)
}
make.variant <- function(external) {
  f <- netVisual_chord_gene; replaced <- c(chord=0L, track=0L, clear=0L)
  rewrite <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name('chordDiagram'))) {x[[1L]] <- as.name('capture.chord'); replaced['chord'] <<- replaced['chord']+1L}
    if (is.call(x) && identical(x[[1L]], as.name('circos.track')) && external) {replaced['track'] <<- replaced['track']+1L; return(quote(plot.labels(lab.cex, grid.col, order.source, order.target)))}
    if (is.call(x) && identical(x[[1L]], as.name('circos.clear')) && external) {replaced['clear'] <<- replaced['clear']+1L; if (replaced['clear']==1L) return(quote({circos.clear(); circlize::circos.par(canvas.xlim=c(-2.1,2.1), canvas.ylim=c(-1.2,1.2))}))}
    if (is.call(x)) return(as.call(lapply(as.list(x), rewrite))); x
  }
  body(f) <- rewrite(body(f)); stopifnot(replaced['chord']==1L, !external || (replaced['track']==1L && replaced['clear']==2L)); f
}
results <- list()
for (current.kind in c('original','external')) {
  f <- make.variant(current.kind == 'external')
  pdf(file.path(out,paste0(current.kind,'.pdf')), width=24, height=16)
  error <- tryCatch({f(chat, signaling=signal, show.legend=FALSE, lab.cex=.7, small.gap=.25); NA_character_}, error=function(e) conditionMessage(e))
  dev.off(); circos.clear()
  results[[current.kind]] <- data.frame(kind=current.kind, error=error, file=file.path(out,paste0(current.kind,'.pdf')))
}
result <- do.call(rbind, results); print(result); write.csv(result,file.path(out,'summary.csv'),row.names=FALSE)
if (any(!is.na(result$error))) stop('At least one prototype plot failed; see summary.csv')
if (all(is.na(result$error))) {
 a <- readRDS(file.path(out,'original-chord-input.rds')); b <- readRDS(file.path(out,'external-chord-input.rds'))
 old.height <- a$args$preAllocateTracks; new.height <- b$args$preAllocateTracks
 a$args$preAllocateTracks <- b$args$preAllocateTracks <- NULL
 stopifnot(identical(a,b)); labels <- readRDS(file.path(out,'labels.rds'))
 print(list(original_track_height=old.height, external_track_height=new.height))
 write.csv(data.frame(links=nrow(a$df), sectors=nrow(labels$positions),
   original_track_height=old.height$track.height, external_track_height=new.height$track.height,
   same_chord_input_except_track_height=TRUE, unique_sector_group_role=TRUE),
   file.path(out,'comparison.csv'), row.names=FALSE)
 for (side in c('left','right')) {
   pos <- labels$positions[labels$positions$side==side,,drop=FALSE]; pos <- pos[order(pos$label.y),,drop=FALSE]
   if (nrow(pos)>1) stopifnot(all(diff(pos$label.y) > (head(pos$label.height,-1)+tail(pos$label.height,-1))/2))
   stopifnot(all(1.22+pos$label.width < 2.1))
 }
 print(c(links=nrow(a$df),sectors=nrow(labels$positions),mapping_rows=nrow(labels$mapping),left=sum(labels$positions$side=='left'),right=sum(labels$positions$side=='right')))
 write.csv(labels$positions,file.path(out,'positions.csv'),row.names=FALSE)
 write.csv(labels$mapping,file.path(out,'mapping.csv'),row.names=FALSE)
}
