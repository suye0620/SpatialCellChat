## Adversarial test: does BiocNeighbors::findNeighbors (KD-tree) miss neighbors
## vs brute-force all-pairs, including boundary distances and duplicate points?
suppressMessages({ library(BiocNeighbors); library(Matrix) })

set.seed(42)
n <- 3000
thr <- 10

## base cloud
xy <- matrix(runif(n * 2, 0, 100), ncol = 2)

## (a) pairs with distance EXACTLY == thr (boundary: original code uses <= thr)
A <- matrix(runif(20 * 2, 0, 90), ncol = 2)
B <- A; B[, 1] <- B[, 1] + thr                    # dist == thr exactly
xy <- rbind(xy, A, B)
n_boundary <- 20

## (b) duplicate coordinate groups (3 copies each x 10 groups)
dup_base <- matrix(runif(10 * 2, 0, 90), ncol = 2)
xy <- rbind(xy, dup_base, dup_base + 1e-12, dup_base + 1e-12)  # near-identical coords
n <- nrow(xy)

## ---- brute force reference (matches current Rfast::Dist logic: 0 < d <= thr) ----
D <- as.matrix(dist(xy))
ref <- lapply(seq_len(n), function(i) which(D[i, ] <= thr & D[i, ] > 0))

## ---- KD-tree radius search ----
kd <- findNeighbors(xy, threshold = thr, get.distance = TRUE)
## exclude self; keep only distance > 0
kd_idx <- lapply(seq_len(n), function(i) {
  idx <- kd$index[[i]]
  d   <- kd$distance[[i]]
  idx[d > 0]          # drop self (d == 0); also drops duplicate-self at d==0
})

## ---- compare, symmetric-normalized ----
norm_pairset <- function(lst) {
  pairs <- unlist(lapply(seq_len(n), function(i) {
    if (length(lst[[i]]) == 0) return(character(0))
    paste(pmin(i, lst[[i]]), pmax(i, lst[[i]]), sep = "-")
  }))
  sort(unique(pairs))
}
pref <- norm_pairset(ref)
pkd  <- norm_pairset(kd_idx)

cat("n =", n, " threshold =", thr, "\n")
cat("ref pairs:", length(pref), " | kd pairs:", length(pkd), "\n")
miss <- setdiff(pref, pkd)
extra <- setdiff(pkd, pref)
cat("MISSING in kd :", length(miss), if (length(miss)) paste(head(miss, 8), collapse = ", ") else "", "\n")
cat("EXTRA in kd   :", length(extra), if (length(extra)) paste(head(extra, 8), collapse = ", ") else "", "\n")

## ---- boundary handling: how many of the exact-distance pairs survive? ----
bd_miss <- sum(grepl(paste0("^", n - 2 * n_boundary + 1, "-"), miss)) +
           sum(grepl(paste0("^", n - n_boundary + 1, "-"), miss))
cat("boundary pairs dropped by kd:", bd_miss, "of", n_boundary, "\n")

## ---- duplicate points: how are d==0 duplicates handled? ----
## brute force keeps all d>0 neighbors; for the dup groups, the second/third copy
## has d ~ 1e-12 > 0 from the first, so should be a neighbor in BOTH.
dup_start <- n - 30 + 1
cat("dup-group first copy neighbor count (ref/kd):",
    length(ref[[dup_start]]), "/", length(kd_idx[[dup_start]]), "\n")
cat("dup-group third copy neighbor count (ref/kd):",
    length(ref[[dup_start + 2]]), "/", length(kd_idx[[dup_start + 2]]), "\n")
