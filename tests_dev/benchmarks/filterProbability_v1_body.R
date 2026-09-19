## 冻结的 filterProbability v1（Wave 1d 原始实现，2026-09-19 从 R/modeling.R L2037-2169 提取）
## 仅对比脚本用；改名避免与 v2 冲突
filterProbability_v1 <- function (
    object,
    nboot = 100,
    seed.use = 666L,
    thresh = 0.05
){
  quantile.prob <- 1-thresh
  if(quantile.prob==0){
    cat(cli.symbol(1),"Do not filter any CCC probability!")
    return(object)
  } else { # quantile.prob<1
    comm.param <- object@misc$.param$communication
    spill.dir <- if (is.null(comm.param)) NULL else comm.param$spill.dir
    if (is.null(object@net$cell$prob) && is.null(spill.dir)) {
      stop(
        cli.symbol(2),
        "Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! "
      )
    } else {
      cat(paste0(cli.symbol(), "Filter out non-significant communication with a probability quantile being ",
                 quantile.prob, " for each L-R pair... \n"))
      # 层来源：内存模式读 net$cell$prob（SparseChatArray 即命名 list）；spill 模式逐层读盘并消费
      if (!is.null(spill.dir)) {
        pair.LR.use <- comm.param$layer.names
        if (is.null(pair.LR.use))
          stop("Spill mode requires `layer.names` in `misc$.param$communication`", call. = FALSE)
        missing.files <- pair.LR.use[!file.exists(file.path(spill.dir, paste0(pair.LR.use, ".rds")))]
        if (length(missing.files))
          stop("Spilled layer files missing for LR pairs: ",
               paste(missing.files, collapse = ", "), call. = FALSE)
        prob.cell_ <- NULL
      } else {
        prob.cell_ <- unclass(object@net$cell$prob)
        pair.LR.use <- dimnames(object@net$cell$prob)[[3]]
      }
      cell.names <- colnames(assay(object, "signaling"))

      d.spatial <- object@images$.distance$d.spatial
      Matrix::diag(d.spatial) <- 1
      adj.contact <- object@images$.distance$adj.contact

      nLR <- length(pair.LR.use)
      nLR1 <- comm.param$nLR1
      nC <- NROW(d.spatial)

      set.seed(seed.use)
      if (comm.param$all.contact.dependent == TRUE) {
        d.spatial <- adj.contact
      }

      # dim(permutation) = nboot x nLR
      permutation <- replicate(nLR, base::sample(x = 1:nC, size = nboot, replace = F))

      prob.cell_ <- my_future_lapply(X = 1:nLR, FUN = function(i) {

        if (i <= nLR1) {
          d_spatial <- d.spatial
        } else {
          d_spatial <- adj.contact
        }
        sample.cells <- permutation[ ,i,drop=T]
        Prob.cell.i <- if (!is.null(spill.dir)) {
          readRDS(file.path(spill.dir, paste0(pair.LR.use[i], ".rds")))
        } else {
          prob.cell_[[i]]
        }

        # 高效分位捷径：基线把各采样源在 pattern 列上的稠密行拼接后取 type-7 分位。
        # 该分位为 0 当且仅当全零前缀长度 n0 满足 n0 >= ceiling(1 + (|S|-1) * quantile.prob)。
        # 零值数下界 = |S| - 层非零值数（O(1) 判定）；不满足时再行级精确判定，
        # 两者均命中则直接返回原层（与基线 nboot.quantile == 0 分支逐位等价），
        # 未命中才构造完整采样向量走基线 quantile 路径。
        pattern.rows <- Matrix::rowSums(d_spatial != 0)
        total.sampled <- sum(pattern.rows[sample.cells])
        cutoff.zero <- FALSE
        if (total.sampled > 0) {
          h.cut <- ceiling(1 + (total.sampled - 1) * quantile.prob)
          if (total.sampled - sum(Prob.cell.i@x != 0) >= h.cut) {
            cutoff.zero <- TRUE
          } else {
            sub.x <- Prob.cell.i[sample.cells, , drop = FALSE]@x
            if (total.sampled - sum(sub.x != 0) >= h.cut) cutoff.zero <- TRUE
          }
        }
        if (cutoff.zero) {
          # nboot.quantile == 0：基线直接返回原层
          return(Prob.cell.i)
        }

        sample.prob.cell.i <- purrr::map(.x = sample.cells,
                                         .f = function(cell.index) {
                                           prob.index <- which(d_spatial[cell.index, ,drop = T] > 0)
                                           nboot.prob.cell.i <- Prob.cell.i[cell.index,prob.index, drop = T] # get a dense vec
                                           return(nboot.prob.cell.i)
                                         }) %>% unlist()
        nboot.quantile <- quantile(sample.prob.cell.i, probs = quantile.prob)
        gc()

        if(nboot.quantile == 0){
          # sparse enough, return directly
          return(Prob.cell.i)
        } else if (nboot.quantile > 0) {
          # filter `Prob.cell.i` to make it sparse enough
          Prob.cell.i <- scMatrixTruncation(Prob.cell.i,cutoff = nboot.quantile,remain.cutoff.v = T,repr = "C")
        }

        return(Prob.cell.i)
      }, simplify = F)
    names(prob.cell_) <- pair.LR.use
    prob.cell <- SparseChatArray(prob.cell_)
    dimnames(prob.cell) <- list(cell.names, cell.names, pair.LR.use)
    object@net$cell$prob <- prob.cell

    if (!is.null(spill.dir)) {
      rds.files <- list.files(spill.dir, pattern = "\\.rds$", full.names = TRUE)
      unlink(rds.files)
      if (length(list.files(spill.dir, all.files = TRUE, no.. = TRUE)) == 0)
        unlink(spill.dir, recursive = TRUE)
      comm.param$spill.dir <- NULL
      comm.param$spill.backend <- NULL
      comm.param$spill.count <- NULL
      comm.param$layer.names <- NULL
      object@misc$.param$communication <- comm.param
    }
    object <- .log_operation(object, "filterProbability", params = list(
      nLR = nLR, nboot = nboot, thresh = thresh,
      spill.consumed = !is.null(spill.dir)))

    cat(cli.symbol(1), "Filtering is done.\n")
    return(object)
    } # whether probability computed
  } # whether to filter out
}
