# 冻结的基线 computeCommunProb 参考实现（2026-09-16 从 R/modeling.R L94-384 逐字转录，
# object 槽位替换为 plain 参数；仅 raw.use 路径）。测试专用 oracle——基线为正确计算基准。
# 修正：computeExpr_complex/_coreceptor 依赖漂移的 my_future_sapply 签名（先在缺陷），
# 此处按其内联语义重写（geometricMean(子矩阵) / ∏(1+e)），数值与基线意图逐位一致。
# 返回 list(prob.cell = sparse3Darray, tmp = list(Lavg, Ravg, prob.cell_), res = 距离缓存)。
baseline_reference_commun_prob <- function(
    data.signaling, LR.use, coordinates, spatial.factors,
    complex_input, cofactor_input, group,
    Kh = 0.5, n = 1, distance.use = TRUE,
    interaction.range = 250, scale.distance = 0.01,
    use.AGAN = TRUE, contact.dependent = TRUE,
    contact.range = 10, contact.dependent.forced = FALSE) {
  data <- data.signaling
  data@x <- data@x / max(data@x)
  data.use <- as.matrix(data)

  if (length(unique(LR.use$annotation)) > 1) {
    LR.use$annotation <- factor(LR.use$annotation,
      levels = c("Secreted Signaling", "ECM-Receptor",
                 "Non-protein Signaling", "Cell-Cell Contact"))
    LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
    LR.use$annotation <- as.character(LR.use$annotation)
  }
  pairLRsig <- LR.use

  ptm <- Sys.time()
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  nC <- ncol(data.use)

  data.spatial <- coordinates
  ratio <- spatial.factors[["ratio"]]
  tol <- spatial.factors[["tol"]]

  res <- computeCellDistance(coordinates = data.spatial, ratio = ratio,
                             interaction.range = interaction.range,
                             contact.range = contact.range, tol = tol)
  d.spatial <- res$d.spatial
  adj.contact <- res$adj.contact

  if (distance.use) {
    d.spatial@x <- d.spatial@x * scale.distance
    d.min <- min(d.spatial@x, na.rm = FALSE)
    if (d.min < 1) stop("scaled distance check failed in reference")
    P.spatial <- createPspatialFrom_dspatial(d.spatial, distance.use = TRUE)
    d.spatial@x <- d.spatial@x / scale.distance
  } else {
    P.spatial <- createPspatialFrom_dspatial(d.spatial, distance.use = FALSE)
  }
  rm(d.spatial); gc()

  all.contact.dependent <- FALSE
  all.diffusible <- FALSE
  if (contact.dependent.forced == TRUE) {
    P.spatial <- P.spatial * adj.contact
    nLR1 <- nLR
    all.contact.dependent <- TRUE
  } else {
    if (contact.dependent == TRUE && length(unique(pairLRsig$annotation)) > 0) {
      if (all(unique(pairLRsig$annotation) == "Cell-Cell Contact")) {
        P.spatial <- P.spatial * adj.contact
        nLR1 <- nLR
        all.contact.dependent <- TRUE
      } else if (all(unique(pairLRsig$annotation) %in%
                     c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
        nLR1 <- nLR
        all.diffusible <- TRUE
      } else {
        nLR1 <- max(which(pairLRsig$annotation %in%
                            c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
      }
    } else {
      nLR1 <- nLR
    }
  }

  # 表达行（基线 computeExpr_LR/_complex 语义内联；复合物 = geometricMean(子矩阵)）
  expr_row <- function(gene) {
    if (gene %in% rownames(data.use)) return(data.use[gene, ])
    sub.cols <- grep("subunit", colnames(complex_input))
    sub <- as.character(unlist(complex_input[gene, sub.cols, drop = FALSE],
                               use.names = FALSE))
    sub <- sub[sub != ""]
    geometricMean(data.use[sub, , drop = FALSE])
  }
  dataLavg <- t(vapply(geneL, expr_row, numeric(nC)))
  dataRavg <- t(vapply(geneR, expr_row, numeric(nC)))
  # 共受体（基线 computeExpr_coreceptor 语义内联：单/多基因统一 ∏(1+e)）
  cof_factor <- function(cof.name) {
    if (is.na(cof.name) || !nzchar(cof.name)) return(rep(1, nC))
    cof.cols <- grep("cofactor", colnames(cofactor_input))
    ind <- cofactor_input[cof.name, cof.cols]
    cof.genes <- as.character(unlist(ind, use.names = FALSE))
    cof.genes <- cof.genes[cof.genes != ""]
    if (!length(cof.genes)) return(rep(1, nC))
    apply(1 + as.matrix(data.use[cof.genes, , drop = FALSE]), 2, prod)
  }
  for (i in seq_len(nLR)) {
    fA <- cof_factor(pairLRsig$co_A_receptor[i])
    fI <- cof_factor(pairLRsig$co_I_receptor[i])
    dataRavg[i, ] <- dataRavg[i, ] * fA / fI
  }

  Prob.cell_ <- lapply(seq_len(nLR), function(i) {
    x_ <- as(matrix(dataLavg[i, ], nrow = 1), "TsparseMatrix")
    y_ <- as(matrix(dataRavg[i, ], nrow = 1), "TsparseMatrix")
    dataLR <- Matrix::crossprod(x_, y_)
    P1 <- HillFunctionFordataLR(dataLR = dataLR, Kh = Kh, n = n)
    P1_Pspatial <- P1 * P.spatial; rm(P1); gc()
    if (!use.AGAN) {
      if (i > nLR1) P1_Pspatial <- P1_Pspatial * adj.contact
      Pnull.cell <- P1_Pspatial
      dimnames(Pnull.cell) <- list(NULL, NULL)
      return(Pnull.cell)
    } else {
      if (i > nLR1) P1_Pspatial <- P1_Pspatial * adj.contact
      if (sum(P1_Pspatial) == 0) {
        Pnull.cell <- P1_Pspatial
        dimnames(Pnull.cell) <- list(NULL, NULL)
        return(Pnull.cell)
      }
      data.agonist <- computeExpr_agonist(data.use, pairLRsig, cofactor_input,
                                          index.agonist = i, Kh = Kh, n = n)
      P_ <- myElementwiseProduct(P1_Pspatial, data.agonist)
      data.antagonist <- computeExpr_antagonist(data.use, pairLRsig, cofactor_input,
                                                index.antagonist = i, Kh = Kh, n = n)
      Pnull.cell <- myElementwiseProduct(P_, data.antagonist)
      rm(P_)
      dimnames(Pnull.cell) <- list(NULL, NULL)
      return(Pnull.cell)
    }
  })

  Prob.cell <- my_as_sparse3Darray(Prob.cell_)
  dimnames(Prob.cell) <- list(colnames(data.use), colnames(data.use), rownames(pairLRsig))
  names(Prob.cell_) <- rownames(pairLRsig)
  Tmp <- list(prob.cell = Prob.cell_, Lavg = dataLavg, Ravg = dataRavg)
  net <- list(prob.cell = Prob.cell, tmp = Tmp)
  list(net = net, res = res,
    run.time = as.numeric(Sys.time() - ptm, units = "secs"),
    parameter = list(Kh = Kh, n = n, nLR = nLR, nLR1 = nLR1,
      scale.distance = scale.distance, use.AGAN = use.AGAN,
      interaction.range = interaction.range,
      all.contact.dependent = all.contact.dependent,
      all.diffusible = all.diffusible))
}
