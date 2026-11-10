# Tests for explicit, display-only arrow shaft scaling.
source('tests_dev/test-grid-visualization.R')
make.plot <- function(scale=1, min.mag=0, ...) suppressWarnings(netVisual_CommunFieldGrid(
  chat,signaling='signal',pattern='outgoing',cellsize=c(1,1),grid.resolution=1,
  arrow.skip=0,min.mag=min.mag,arrow.scale=scale,...))
input.before<-chat
p.default<-suppressWarnings(netVisual_CommunFieldGrid(chat,signaling='signal',pattern='outgoing',cellsize=c(1,1),grid.resolution=1,arrow.skip=0))
p.one<-make.plot(1);p.four<-make.plot(4)
a<-ggplot_build(p.default);b<-ggplot_build(p.one);c<-ggplot_build(p.four)
stopifnot(identical(chat,input.before),identical(a$data,b$data),
  is.null(p.default$labels$caption),is.null(p.one$labels$caption),
  identical(p.four$labels$caption,'Arrow length display scale: 4x'),
  identical(p.default$layers[[2L]]$data,p.four$layers[[2L]]$data),
  identical(a$data[[1L]],c$data[[1L]]))
common<-setdiff(intersect(names(a$data[[2L]]),names(c$data[[2L]])),'mag')
stopifnot(isTRUE(all.equal(a$data[[2L]][,common],c$data[[2L]][,common],tolerance=0,check.attributes=TRUE)),
  isTRUE(all.equal(c$data[[2L]]$mag,a$data[[2L]]$mag*4,tolerance=0)))
# Verify the grid grob actually draws fourfold shaft lengths while arrowheads retain their prior settings.
arrow.grob<-function(plot) {
  g<-ggplot2::ggplotGrob(plot);panel<-g$grobs[[which(g$layout$name=='panel')]]
  found<-Filter(function(x)inherits(x,'arrow2'),panel$children)
  if(length(found)!=1L)stop('Expected one metR arrow grob')
  found[[1L]]
}
grob.one<-arrow.grob(p.one);grob.four<-arrow.grob(p.four)
stopifnot(isTRUE(all.equal(as.numeric(grob.four$length),4*as.numeric(grob.one$length),tolerance=0)),
  identical(grob.one$arrow,grob.four$arrow),identical(grob.one$angle,grob.four$angle))
# Threshold selection is based on original, not display-scaled, magnitudes.
threshold<-max(a$data[[2L]]$mag,na.rm=TRUE)*.5
threshold.one<-ggplot_build(make.plot(1,min.mag=threshold))$data[[2L]]
threshold.four<-ggplot_build(make.plot(4,min.mag=threshold))$data[[2L]]
stopifnot(nrow(threshold.one)>0,nrow(threshold.one)<nrow(a$data[[2L]]),
  isTRUE(all.equal(threshold.four$mag,4*threshold.one$mag,tolerance=0)),
  identical(threshold.one$x,threshold.four$x),identical(threshold.one$y,threshold.four$y))
custom<-make.plot(4,pivot=0)
stopifnot(identical(custom$layers[[2L]]$geom_params$pivot,0))
for(value in list(0,-1,Inf,NA_real_,numeric(),c(1,2),'four',TRUE,1+1i)) {
  error<-tryCatch({make.plot(value);''},error=conditionMessage)
  if(!grepl('arrow.scale must be one positive finite number',error,fixed=TRUE)) stop('Invalid value not rejected: ',paste(value,collapse=','),'; ',error)
}
zero.scaled<-suppressWarnings(netVisual_CommunFieldGrid(zero.chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,arrow.skip=0,arrow.scale=4))
stopifnot(all(ggplot_build(zero.scaled)$data[[2L]]$alpha==0))
pdf<-tempfile(fileext='.pdf');rds<-tempfile(fileext='.rds');ggsave(pdf,p.four,width=6,height=5);saveRDS(p.four,rds)
stopifnot(file.info(pdf)$size>0,isTRUE(all.equal(ggplot_build(readRDS(rds))$data[[2L]],c$data[[2L]],tolerance=0)))
unlink(c(pdf,rds));cat('Arrow display-scale checks passed.\n')
