# Archived: superseded expression / aggregation helpers
#
# Archived on 2026-10-26 during the 11-slot migration. These functions are no
# longer part of the package: R only compiles top-level .R files under R/, so
# nothing in this directory is sourced by R CMD INSTALL. They are kept here for
# reference and for the frozen baseline comparison scripts in tests_dev/.
#
# Why they were retired:
#   - the per-pair expression gating moved into the Rcpp kernels
#     (cpp_prob_layer, cpp_group_avg_obs, cpp_group_avg_perm);
#   - the 3-D array container moved from sparse3Darray to SparseChatArray;
#   - cofactor / co-receptor / agonist-antagonist modelling now goes through
#     .sc_expr_row(), .sc_complex_subunits() and .sc_cofactor_factor().
#
# Original locations (line numbers at archive time):
#   R/modeling.R   483  HillFunctionFordataLR
#   R/modeling.R   505  elementwiseProductForSparseMatandDenseMat
#   R/modeling.R   536  myElementwiseProduct
#   R/modeling.R  1378  computeAvgCommunProb_LR_Sum
#   R/modeling.R  1449  computeAvgCommunProb_LR_Avg
#   R/modeling.R  2046  modifySparse3Darray
#   R/modeling.R  2565  .computeExprGroup_complex
#   R/modeling.R  2631  computeExpr_coreceptor
#   R/modeling.R  2678  .computeExprGroup_coreceptor
#   R/modeling.R  2737  computeExprGroup_agonist
#   R/modeling.R  2769  computeExprGroup_antagonist
#   R/modeling.R  2803  computeExpr_agonist
#   R/modeling.R  2838  computeExpr_antagonist
#   R/utilities.R  361  ==.sparse3Darray
#   R/utilities.R  491  my_as_sparse3Darray
#
# Restore: move the block back to the file named above and re-add its
# export() line in NAMESPACE (the 12 exported ones are listed in the agent
# note .agents/notes/implemented/2026-10-26-archive-expression-toolchain.md).
#
# KNOWN CONSEQUENCE: netVisual_CommunField() (R/visualization.R) still calls
# my_as_sparse3Darray() on the retired net$field slot. That path was already
# unreachable before this archive; the function is queued for the next
# archiving round together with netVisual_CommunFieldZoomIn().


# ==========================================================================
# Hill 函数与 element-wise 乘积（旧 computeCommunProb 的内循环实现）
# ==========================================================================

#' HillFunction for dataLR object
#' @description
#' A function defined for \link{computeCommunProb}. A Hill function was used to model the interactions
#' between L and R with an EC50 parameter Kh whose default value was set to be 0.5 as the input data has a normalized range from 0 to 1.
#' Previous hill function is `dataLR^n/(Kh^n + dataLR^n)`, notice dataLR is a sparse matrix, so we design the function to do element-wise hill function upon each non-zero value in the `dataLR` matrix
#'
#' @param dataLR A dgCMatrix object in `Matrix` Package.
#' @param Kh A hyper parameter in hill function, passed from \link{computeCommunProb}
#' @param n A hyper parameter in hill function, passed from \link{computeCommunProb}
#'
#' @return a normalized dataLR
#' @export
#'
#' @examples
HillFunctionFordataLR <- function(dataLR,Kh,n){
  # dataLR is a dgCMatrix, we map the Hill function on each non-zero in the Matrix
  x_dataLR_ <- purrr::map_dbl(dataLR@x,function(x){x^n/(Kh^n + x^n)})
  # update the values
  dataLR@x <- x_dataLR_
  return(dataLR)
}

#' element-wise product for a sparse matrix and a dense matrix
#'
#' @description
#' A function defined for \link{computeCommunProb}. Define an element-wise product operation for a sparse matrix and a dense matrix,
#' to speed up the product operation
#'
#' @param SparseMat A dgCMatrix object in `Matrix` Package.
#' @param DenseMat A dense matrix shares the same dimensions with `SparseMat`
#'
#' @return element-wise product result. A dgCMatrix object in `Matrix` Package
#' @export
#'
#' @examples
elementwiseProductForSparseMatandDenseMat <- function(SparseMat,DenseMat){
  # SparseMat must be a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "TsparseMatrix")
  nonZeroNum <- length(SparseMat@i)

  xFromDenseMat <- vector('double',length = nonZeroNum)
  for (k in 1:nonZeroNum) {
    # row indices of non-zero entries in 0-base, so add 1
    # col indices of non-zero entries in 0-base, so add 1
    xFromDenseMat[[k]] <- DenseMat[(SparseMat@i[[k]]+1), (SparseMat@j[[k]]+1)]
  }
  # update the x in SparseMat
  SparseMat@x <- SparseMat@x * xFromDenseMat
  # still return a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "CsparseMatrix")
  return(SparseMat)
}

#' element-wise product for a sparse matrix and a dense vector's crossprod
#'
#' @description
#' A function defined for \link{computeCommunProb}. Define an element-wise product operation for a sparse matrix and a dense vector's crossprod,
#' to speed up the product operation
#'
#' @param SparseMat A dgCMatrix object in `Matrix` Package.
#' @param DenseVec A numeric vector shares the same dimensions with `SparseMat`
#'
#' @return element-wise product result. A dgCMatrix object in `Matrix` Package
#' @export
myElementwiseProduct <- function(SparseMat,DenseVec){
  # SparseMat must be a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "TsparseMatrix")
  nonZeroNum <- length(SparseMat@i)

  xFromCrossProduct <- vector('double',length = nonZeroNum)
  for (k in 1:nonZeroNum) {
    # row indices of non-zero entries in 0-base, so add 1
    # col indices of non-zero entries in 0-base, so add 1
    # DenseVec[[(SparseMat@i[[k]]+1)]] * DenseVec[[(SparseMat@j[[k]]+1)]] = crossprod[i,j]
    xFromCrossProduct[[k]] <- DenseVec[[(SparseMat@i[[k]]+1)]] * DenseVec[[(SparseMat@j[[k]]+1)]]
  }
  # update the x in SparseMat
  SparseMat@x <- SparseMat@x * xFromCrossProduct
  # still return a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "CsparseMatrix")
  return(SparseMat)
}



# ==========================================================================
# 旧 group-average 聚合（被 computeAvgCommunProb 的 Rcpp 核取代）
# ==========================================================================

#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
#' Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_Sum <- function (
    prob,
    group,
    dataLR = NULL,
    min.percent = 0.1,
    min.cells.sr = 5
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5

  # group is vector/factor
  if (!is.factor(group)) {
    group <- factor(group)
  }
  cell.type.mat <- model.matrix(~group-1)


  # min.percent
  dataLR_temp <- 1 * (dataLR > 0)
  dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
  dataLR_percent <- 1 * (signif(dataLR_temp[, -1], 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))

  if(sum(Prob_percent)==0){
    Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
  } else {

    Prob.avg <- Matrix::crossprod(x = cell.type.mat,
                                  y = prob %*% cell.type.mat)

    # min.cells.sr
    prob@x <- rep.int(1, times = length(prob@x))
    sender.counts <- Matrix::rowSums(prob)
    receptor.counts <- Matrix::colSums(prob)
    cells.sr <- cbind(sender.counts, receptor.counts)
    cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
    cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))

    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }

  # return
  Prob.avg <- as.matrix(Prob.avg)
  return(Prob.avg)
}

#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param min.percent Numeric from 0 to 1. The expressed-ligand/receptor proportion is rounded to one significant digit numerically and must meet this threshold for group-level signaling
#' Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_Avg <- function (
    prob,
    group,
    dataLR = NULL,
    min.percent = 0.1,
    min.cells.sr = 5
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5

  # group is vector/factor
  if (!is.factor(group)) {
    group <- factor(group)
  }
  cell.type.mat <- model.matrix(~group-1)


  # min.percent
  dataLR_temp <- 1 * (dataLR > 0)
  dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
  dataLR_percent <- 1 * (signif(dataLR_temp[, -1], 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))

  if(sum(Prob_percent)==0){
    Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
  } else {
    Prob.avg <- Matrix::crossprod(x = cell.type.mat,
                                  y = prob %*% cell.type.mat)

    # min.cells.sr
    prob@x <- rep.int(1, times = length(prob@x))
    Prob.scale.factor <- Matrix::crossprod(x = cell.type.mat,
                                           y = prob %*% cell.type.mat)
    Prob.avg <- Prob.avg/Prob.scale.factor
    Prob.avg[is.nan(Prob.avg)] <- 0

    sender.counts <- Matrix::rowSums(prob)
    receptor.counts <- Matrix::colSums(prob)
    cells.sr <- cbind(sender.counts, receptor.counts)
    cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
    cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))

    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }
  Prob.avg <- as.matrix(Prob.avg)

  return(Prob.avg)
}



# ==========================================================================
# 旧 sparse3Darray 桥接（被 SparseChatArray 取代）
# ==========================================================================

#' @title modifySparse3Darray
#'
#' @param array3d sparse3Darray. See details in \code{\link[spatstat.sparse]{sparse3Darray}}
#' @param do.binary Boolean. Use 1 to replace the non-zero values in the 3Darray.
#' @param cutoff Numeric. Cut off the sparse3Darray's values to make
#' array' values >= cutoff.
#' @param remain.cutoff.v Boolean. Use `>=` or `>`.
#'
#' @return sparse3Darray.
#' @export
modifySparse3Darray <- function(
    array3d,
    do.binary=T,
    cutoff=NULL,
    remain.cutoff.v=T
){
  if (!inherits(array3d, "sparse3Darray")) {
    stop("Please use sparse3Darray as input!")
  }

  if((!is.null(cutoff)) & is.numeric(cutoff)){
    if(length(array3d$x) > 0){

      if(remain.cutoff.v){
        remain.idx <- which(array3d$x >= cutoff)
      } else {
        remain.idx <- which(array3d$x > cutoff)
      }

      # update
      array3d$i <- array3d$i[remain.idx]
      array3d$j <- array3d$j[remain.idx]
      array3d$k <- array3d$k[remain.idx]
      array3d$x <- array3d$x[remain.idx]
    }
  }

  if(do.binary){
    array3d$x <- rep.int(1,times = length(array3d$x))
  }

  return(array3d)
}

#' A faster realization in comparison to \code{\link[spatstat.sparse]{as.sparse3Darray}}
#'
#' @param x Data in another format (see Details in \code{\link[spatstat.sparse]{as.sparse3Darray}}.
#' @param ... Ignored.
#' @param strict Boolean. Logical value specifying whether to enforce the rule that each entry in i,j,k,x refers to a different cell. If strict=TRUE, entries which refer to the same cell in the array will be reduced to a single entry by summing the x values. Default is strict=FALSE.
#' @param nonzero Boolean. Logical value specifying whether to remove any entries of x which equal zero.
#' @importFrom spatstat.sparse sparse3Darray as.sparse3Darray
#' @return Sparse three-dimensional array (object of class "sparse3Darray").
#' @export
#'
my_as_sparse3Darray <- function (x, ... ,strict = FALSE, nonzero = FALSE)
{
  if (inherits(x, "sparse3Darray")) {
    y <- x
  }
  else if (inherits(x, c("matrix", "sparseMatrix"))) {
    z <- as(x, Class = "TsparseMatrix")
    dn <- dimnames(x)
    dn <- if (is.null(dn))
      NULL
    else c(dn, list(NULL))
    one <- if (length(z@i) > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = z@i + 1L, j = z@j + 1L, k = one,
                       x = z@x, dims = c(dim(x), 1L), dimnames = dn,strict = strict, nonzero = nonzero)
  }
  else if (is.array(x)) {
    stopifnot(length(dim(x)) == 3)
    dimx <- dim(x)
    if (prod(dimx) == 0) {
      y <- sparse3Darray(, dims = dimx, dimnames = dimnames(x))
    }
    else {
      ijk <- which(x != spatstat.utils::RelevantZero(x), arr.ind = TRUE)
      ijk <- cbind(as.data.frame(ijk), x[ijk])
      y <- sparse3Darray(i = ijk[, 1L], j = ijk[, 2L],
                         k = ijk[, 3L], x = ijk[, 4L], dims = dimx, dimnames = dimnames(x),strict = strict, nonzero = nonzero)
    }
  }
  else if (inherits(x, "sparseVector")) {
    one <- if (length(x@i) > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = x@i, j = one, k = one, x = x@x,
                       dims = c(x@length, 1L, 1L),strict = strict, nonzero = nonzero)
  }
  else if (is.null(dim(x)) && is.atomic(x)) {
    n <- length(x)
    dn <- names(x)
    if (!is.null(dn))
      dn <- list(dn, NULL, NULL)
    one <- if (n > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = seq_len(n), j = one, k = one,
                       x = x, dims = c(n, 1L, 1L), dimnames = dn,strict = strict, nonzero = nonzero)
  }
  else if (is.list(x) && length(x) > 0) {
    n <- length(x)
    if (all(sapply(x, is.matrix))) {
      z <- Reduce(abind, x)
      y <- as.sparse3Darray(z)
    }
    else if (all(sapply(x, inherits, what = "sparseMatrix"))) {
      dimlist <- unique(lapply(x, dim))
      if (length(dimlist) > 1)
        stop("Dimensions of matrices do not match")
      dimx <- c(dimlist[[1L]], n)
      dnlist <- lapply(x, dimnames)
      isnul <- sapply(dnlist, is.null)
      dnlist <- unique(dnlist[!isnul])
      if (length(dnlist) > 1)
        stop("Dimnames of matrices do not match")
      dn <- if (length(dnlist) == 0)
        NULL
      else c(dnlist[[1L]], list(NULL))

      progressr::handlers(global = T)
      progressr::handlers(progressr::handler_progress(format = ":percent [:bar] :eta :message"))
      X <- seq_len(n)
      p <- progressr::progressor(along = X)
      df <- purrr::map_dfr(
        .x = X,
        .f = function(k) {
          mk <- as(x[[k]], "TsparseMatrix")
          kvalue <- if (length(mk@i) > 0) k else integer(0)
          dfk <- data.frame(
            i = mk@i + 1L,
            j = mk@j + 1L,
            k = kvalue,
            x = mk@x
          )
          p("as sparse3Darray...")
          return(dfk)
        }
      )

      y <- sparse3Darray(i = df$i, j = df$j, k = df$k,
                         x = df$x, dims = dimx, dimnames = dn,strict = strict, nonzero = nonzero)
    }
    else {
      warning("I don't know how to convert a list to a sparse array")
      return(NULL)
    }
  }
  else {
    warning("I don't know how to convert x to a sparse array")
    return(NULL)
  }
  return(y)
}

#' @title Fast `==` operation for `sparse3Darray` Class
#'
#' @param arr1 A sparse3Darray
#' @param arr2 A sparse3Darray
#'
#' @return Boolean.
#' @export
`==.sparse3Darray` <- function(arr1, arr2) {
  if (!inherits(arr1, "sparse3Darray") || !inherits(arr2, "sparse3Darray")) {
    stop("Both arguments must be of class 'sparse3Darray'")
  }
  i.equal = all(arr1$i==arr2$i)
  j.equal = all(arr1$j==arr2$j)
  k.equal = all(arr1$k==arr2$k)
  x.equal = all(arr1$x==arr2$x)
  dim.equal = all(arr1$dim==arr2$dim)

  return(all(c(i.equal,j.equal,k.equal,x.equal,dim.equal)))
}



# ==========================================================================
# 旧表达 / 共受体 / 激动-拮抗剂建模（被 .sc_expr_row 与 .sc_complex_subunits 取代）
# ==========================================================================

# Compute the average expression of complex per cell group using geometric mean
# @param complex_input the complex_input from CellChatDB
# @param data.use data matrix (rows are genes and columns are cells)
# @param complex the names of complex
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom dplyr select starts_with
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_complex <- function(complex_input, data.use, complex, group, FunMean) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = sapply,
    no = future.apply::future_sapply
  )
  data.complex = my.sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      RsubunitsV <- intersect(RsubunitsV, rownames(data.use))
      if (length(RsubunitsV) > 1) {
        data.avg <- aggregate(t(data.use[RsubunitsV,]), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else if (length(RsubunitsV) == 1) {
        data.avg <- aggregate(matrix(data.use[RsubunitsV,], ncol = 1), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else {
        data.avg = matrix(0, nrow = 1, ncol = length(unique(group)))
      }
      return(geometricMean(data.avg))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

#' Modeling the effect of coreceptor on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig a data frame giving ligand-receptor interactions
#' @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
#' @return
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I")) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {

    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my_future_sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) == 1) {
          return(1 + data.use[coreceptor.indV, ])
        } else if (length(coreceptor.indV) > 1) {
          return(apply(1 + data.use[coreceptor.indV, ], 2, prod))
        } else {
          return(matrix(1, nrow = 1, ncol = ncol(data.use)))
        }
      },
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
  }
  return(data.coreceptor)
}

# Modeling the effect of coreceptor on the ligand-receptor interaction
#
# @param data.use data matrix
# @param cofactor_input the cofactor_input from CellChatDB
# @param pairLRsig a data frame giving ligand-receptor interactions
# @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I"), group, FunMean) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = pbapply::pbsapply,
      no = future.apply::future_sapply
    )
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my.sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) > 1) {
          data.avg <- aggregate(t(data.use[coreceptor.indV,]), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(apply(1 + data.avg, 2, prod))
          # return(1 + apply(data.avg, 2, mean))
        } else if (length(coreceptor.indV) == 1) {
          data.avg <- aggregate(matrix(data.use[coreceptor.indV,], ncol = 1), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(1 + data.avg)
        } else {
          return(matrix(1, nrow = 1, ncol = length(unique(group))))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
  }

  return(data.coreceptor)
}

#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
#' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_agonist <- function(data.use, pairLRsig, cofactor_input, group, index.agonist, Kh, FunMean, n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
#' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_antagonist <- function(data.use, pairLRsig, cofactor_input, group, index.antagonist, Kh, FunMean, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.antagonist)
}

#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
# #' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_agonist <- function(data.use, pairLRsig, cofactor_input, index.agonist, Kh,  n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.agonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
# #' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_antagonist <- function(data.use, pairLRsig, cofactor_input, index.antagonist, Kh, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.antagonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.antagonist)
}
