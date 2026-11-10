# Diagnose display sampling and streamline integration on the full Slide-seq fixture.
setwd(local({
  args <- commandArgs(FALSE)
  file <- sub('^--file=', '', grep('^--file=', args, value = TRUE))
  if (length(file)) dirname(dirname(normalizePath(file))) else getwd()
}))
Sys.setenv(RENV_PATHS_LIBRARY = 'renv/library')
if (!nzchar(Sys.getenv('RENV_PROJECT'))) {
  if (requireNamespace('renv', quietly = TRUE)) renv::load(getwd()) else source('renv/activate.R')
}
suppressPackageStartupMessages({
  library(methods); library(Matrix); library(cli); library(dplyr); library(tidyr)
  library(purrr); library(sf); library(ggplot2); library(colorspace); library(spatstat.sparse)
})
source('R/SpatialCellChat_class.R'); source('R/database.R'); source('R/utilities.R')
source('R/spatial.R'); source('R/modeling.R'); source('R/analysis.R'); source('R/visualization.R')
source('tests_dev/fixtures/make_slideseq_window.R')
chat <- make_slideseq_chat()
signal <- dimnames(chat@netP$cell$prob)[[3L]][1L]
chat <- computeCommunField(chat, signaling.name = signal)
field <- as.matrix(chat@netP$cell$field$outgoing[[signal]])
out.dir <- 'tests_dev/figures/P3/slideseq/display-diagnostics-2026-11-11'
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
warnings <- list()
capture.warnings <- function(name, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- data.frame(case = name, warning = conditionMessage(w))
    invokeRestart('muffleWarning')
  })
}
grid <- capture.warnings('grid-construction', netVisual_CommunFieldGrid(chat, signaling = signal, pattern = 'outgoing'))
flow <- capture.warnings('flow-construction', netVisual_CommunFlow(chat, signaling = signal, pattern = 'outgoing'))
uv <- grid$layers[[2L]]$data
stopifnot(nrow(field) == nrow(chat@meta), identical(rownames(field), rownames(chat@images$coordinates)))
# Sampling and opacity changes are restricted to the already computed plot layer.
copy.layers <- function(plot) {
  plot$layers <- lapply(plot$layers, function(layer) ggplot2::ggproto(NULL, layer))
  plot
}
settings <- list(grid_default = copy.layers(grid), grid_no_skip = copy.layers(grid),
                 flow_default = copy.layers(flow), flow_res8 = copy.layers(flow))
settings$grid_no_skip$layers[[2L]]$stat_params$skip.x <- 0L
settings$grid_no_skip$layers[[2L]]$stat_params$skip.y <- 0L
settings$grid_no_skip$layers[[1L]]$aes_params$size <- 0.22
settings$grid_no_skip$layers[[1L]]$aes_params$alpha <- 0.4
settings$flow_res8$layers[[2L]]$stat_params$res <- 8
settings$flow_res8$layers[[1L]]$aes_params$size <- 0.22
settings$flow_res8$layers[[1L]]$aes_params$alpha <- 0.4
stopifnot(settings$grid_default$layers[[2L]]$stat_params$skip.x == 4L,
          settings$flow_default$layers[[2L]]$stat_params$res == 0.5,
          settings$grid_default$layers[[1L]]$aes_params$size == 2,
          settings$flow_default$layers[[1L]]$aes_params$alpha == 0.55)
stats <- list()
for (name in names(settings)) {
  p <- settings[[name]]
  built <- capture.warnings(paste0(name, '-build'), ggplot2::ggplot_build(p))
  drawn <- built$data[[2L]]
  params <- p$layers[[2L]]$stat_params
  file <- file.path(out.dir, paste0(name, '.pdf'))
  capture.warnings(paste0(name, '-pdf'), ggplot2::ggsave(file, p, width = 7, height = 5.5))
  if (!file.exists(file) || file.info(file)$size < 1000L) stop('missing or empty PDF: ', file)
  is.flow <- startsWith(name, 'flow')
  group.sizes <- if (is.flow && nrow(drawn)) table(drawn$group) else integer()
  group.sizes <- group.sizes[group.sizes > 0L]
  integration <- if (is.flow) {
    input <- p$layers[[2L]]$data
    input <- data.frame(x = input$lon, y = -input$lat, dx = input$u, dy = -input$v, group = 1L)
    suppressWarnings(metR::StatStreamline$setup_params(input, params))
  } else NULL
  stats[[length(stats) + 1L]] <- data.frame(
    case = name, dataset = 'slideseq', signal = signal, cells = nrow(field),
    nonzero_field_cells = sum(rowSums(field^2) > 0),
    interpolated_points = nrow(uv), nonzero_interpolated_points = sum(uv$u^2 + uv$v^2 > 0, na.rm = TRUE),
    arrow_skip = if (is.flow) NA_integer_ else params$skip.x,
    flow_res = if (is.flow) params$res else NA_real_,
    flow_length = if (is.flow) params$L else NA_real_,
    integration_steps_each_direction = if (is.flow) integration$S else NA_integer_,
    integration_step_dt = if (is.flow) integration$dt else NA_real_,
    drawn_rows = nrow(drawn),
    visible_positive_arrows = if (is.flow) NA_integer_ else sum(drawn$mag > 0 & drawn$alpha > 0, na.rm = TRUE),
    zero_transparent_arrows = if (is.flow) NA_integer_ else sum(drawn$mag == 0 & drawn$alpha == 0, na.rm = TRUE),
    median_positive_arrow_cm = if (is.flow || !any(drawn$mag > 0, na.rm = TRUE)) NA_real_ else median(drawn$mag[drawn$mag > 0], na.rm = TRUE),
    max_positive_arrow_cm = if (is.flow || !any(drawn$mag > 0, na.rm = TRUE)) NA_real_ else max(drawn$mag[drawn$mag > 0], na.rm = TRUE),
    streamline_groups = if (is.flow) length(group.sizes) else NA_integer_,
    streamline_groups_at_least_3_points = if (is.flow) sum(group.sizes >= 3) else NA_integer_,
    point_size = p$layers[[1L]]$aes_params$size, point_alpha = p$layers[[1L]]$aes_params$alpha,
    file = file)
}
summary <- do.call(rbind, stats)
write.csv(summary, file.path(out.dir, 'summary.csv'), row.names = FALSE)
warning.table <- if (length(warnings)) unique(do.call(rbind, warnings)) else data.frame(case = character(), warning = character())
write.csv(warning.table, file.path(out.dir, 'warnings.csv'), row.names = FALSE)
print(summary); print(warning.table)
