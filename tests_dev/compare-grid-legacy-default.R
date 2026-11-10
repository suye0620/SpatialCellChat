# Compare current defaults with the checked-in pre-option implementation. Requires Git history.
source('tests_dev/test-grid-visualization.R')
old.file<-tempfile(fileext='.R');status<-system2('git',c('show','66a16b6:R/visualization.R'),stdout=old.file)
stopifnot(identical(status,0L));old.env<-new.env(parent=.GlobalEnv);sys.source(old.file,envir=old.env)
old.plot<-suppressWarnings(old.env$netVisual_CommunFieldGrid(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,arrow.skip=0))
new.plot<-suppressWarnings(netVisual_CommunFieldGrid(chat,signaling='signal',pattern='outgoing',cellsize=1,grid.resolution=1,arrow.skip=0))
old.built<-ggplot2::ggplot_build(old.plot);new.built<-ggplot2::ggplot_build(new.plot)
stopifnot(identical(old.plot$layers[[2L]]$data,new.plot$layers[[2L]]$data),identical(old.built$data,new.built$data),identical(old.plot$labels,new.plot$labels))
unlink(old.file);cat('Unmodified pre-change function and current default: identical layer inputs, built layers, and labels.\n')
