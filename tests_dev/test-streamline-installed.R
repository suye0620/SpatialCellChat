# Test the installed namespace and reloading its returned plot in a fresh R process.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('Pass the temporary package library')
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
suppressPackageStartupMessages(library(SpatialCellChat))
coords <- as.matrix(expand.grid(x=0:4,y=0:4)); cells <- paste0('c',seq_len(nrow(coords))); rownames(coords) <- cells
expr <- Matrix::Matrix(matrix(1,nrow=1,ncol=length(cells),dimnames=list('G1',cells)),sparse=TRUE)
meta <- data.frame(label=factor(rep('A',length(cells))),row.names=cells)
chat <- createSpatialCellChat(expr,meta=meta,group.by='label',datatype='spatial',coordinates=coords,spatial.factors=list(ratio=1,tol=.5))
components <- matrix(rep(c(1,1),each=length(cells)),ncol=2,dimnames=list(cells,c('x','y')))
chat@netP$cell$field$outgoing <- SpatialCellChat:::SparseChatArray(list(signal=as(components,'dgCMatrix')),dimnames=list(cells,c('x','y'),'signal'))
p <- netVisual_CommunFlow(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0)
stopifnot(nrow(ggplot2::ggplot_build(p)$data[[2L]])>0)
file <- tempfile(fileext='.rds'); pdf.file <- tempfile(fileext='.pdf'); child <- tempfile(fileext='.R')
saveRDS(p,file)
writeLines(c('.libPaths(c(commandArgs(TRUE)[1L], .libPaths()))',
  'p <- readRDS(commandArgs(TRUE)[2L])',
  'if (!nrow(ggplot2::ggplot_build(p)$data[[2L]])) stop("No streamlines after reload")',
  'ggplot2::ggsave(commandArgs(TRUE)[3L],p,width=6,height=5)'),child)
status <- system2(file.path(R.home('bin'),'Rscript.exe'),c(shQuote(child),shQuote(normalizePath(args[[1L]])),shQuote(file),shQuote(pdf.file)))
stopifnot(identical(status,0L),file.info(pdf.file)$size>0)
unlink(c(file,pdf.file,child))
cat('Installed namespace and fresh-process serialization checks passed.\n')
