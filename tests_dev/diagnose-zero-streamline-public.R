# Reproduce all-zero flow behavior through the public entry point.
setwd(local({a<-commandArgs(FALSE); f<-sub('^--file=','',grep('^--file=',a,value=TRUE));if(length(f))dirname(dirname(normalizePath(f)))else getwd()}))
Sys.setenv(RENV_PATHS_LIBRARY='renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {if (requireNamespace('renv',quietly=TRUE)) renv::load(getwd()) else source('renv/activate.R')}
suppressPackageStartupMessages({library(methods);library(Matrix);library(cli);library(dplyr);library(tidyr);library(purrr);library(sf);library(ggplot2);library(colorspace);library(spatstat.sparse)})
source('R/SpatialCellChat_class.R');source('R/database.R');source('R/utilities.R');source('R/spatial.R');source('R/modeling.R');source('R/analysis.R');source('R/visualization.R')
coords<-as.matrix(expand.grid(x=0:4,y=0:4));cells<-paste0('c',seq_len(nrow(coords)));rownames(coords)<-cells
expr<-Matrix::Matrix(matrix(1,nrow=1,ncol=length(cells),dimnames=list('G1',cells)),sparse=TRUE)
meta<-data.frame(label=factor(rep('A',length(cells))),row.names=cells)
chat<-createSpatialCellChat(expr,meta=meta,group.by='label',datatype='spatial',coordinates=coords,spatial.factors=list(ratio=1,tol=.5))
components<-Matrix::Matrix(matrix(0,nrow=length(cells),ncol=2,dimnames=list(cells,c('x','y'))),sparse=TRUE)
chat@netP$cell$field$outgoing<-SparseChatArray(list(signal=as(components,'dgCMatrix')),dimnames=list(cells,c('x','y'),'signal'))
warnings<-character();p<-suppressWarnings(netVisual_CommunFlow(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,L=2,min.L=0,res=8,start=data.frame(x=2,y=-2),jitter=0,stat=ggplot2::ggproto(NULL,metR::StatStreamline)))
value<-tryCatch(withCallingHandlers(ggplot2::ggplot_build(p)$data[[2L]],warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
cat('interpolated rows:',nrow(p$layers[[2L]]$data),' nonzero:',sum(p$layers[[2L]]$data$u!=0|p$layers[[2L]]$data$v!=0,na.rm=TRUE),' result rows:',if(inherits(value,'error')) NA_integer_ else nrow(value),'\n')
if(inherits(value,'error'))cat('error:',conditionMessage(value),'\n')
if(packageVersion('metR') == package_version('0.19.0')) stopifnot(inherits(value,'error'),grepl('missing value where TRUE/FALSE needed',conditionMessage(value),fixed=TRUE))
cat(paste(warnings,collapse='\n'),'\n')

# Prototype a zero-row drawing layer, leaving the original plot and stored field untouched.
q<-p;q$layers[[2L]]<-ggplot2::ggproto(NULL,p$layers[[2L]]);q$layers[[2L]]$data<-q$layers[[2L]]$data[FALSE,,drop=FALSE]
fixed<-ggplot2::ggplot_build(q)$data[[2L]]
cat('empty drawing layer rows:',nrow(fixed),' background rows:',nrow(ggplot2::ggplot_build(q)$data[[1L]]),'\n')
stopifnot(nrow(fixed)==0L,nrow(ggplot2::ggplot_build(q)$data[[1L]])==length(cells))
out<-tempfile(fileext='.pdf');ggplot2::ggsave(out,q,width=6,height=5);stopifnot(file.info(out)$size>0);unlink(out)
