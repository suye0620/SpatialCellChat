#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

typedef Eigen::Triplet<double> T;
// Adapted from swne (https://github.com/yanwu2014/swne)
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}

// Sum N identically-dimensioned dgCMatrix into a single dgCMatrix.
// Column-by-column accumulator: val[n]+mark[n] per column avoids intermediate expansion.
// [[Rcpp::export]]
S4 cpp_sum_layers(List matrices) {
  int n = matrices.size();
  if (n == 0) stop("matrices is empty");

  // Extract slots from every matrix up front
  int nr = 0, nc = 0;
  std::vector<IntegerVector> all_p(n), all_i(n);
  std::vector<NumericVector> all_x(n);
  for (int k = 0; k < n; k++) {
    S4 M(matrices[k]);
    IntegerVector d = M.slot("Dim");
    if (k == 0) { nr = d[0]; nc = d[1]; }
    else if (d[0] != nr || d[1] != nc)
      stop("all matrices must have the same dimensions");
    all_p[k] = M.slot("p");
    all_i[k] = M.slot("i");
    all_x[k] = M.slot("x");
  }

  // Column-local accumulators
  std::vector<double> val(nr, 0.0);
  std::vector<int>    mark(nr, -1);   // -1 = untouched; otherwise column id in current pass
  std::vector<int>    rows;           // row indices touched in current column
  rows.reserve(nr);

  // ---- Pass 1: count non-zeros per column (deduplicated) ----
  IntegerVector p_out(nc + 1);
  p_out[0] = 0;
  for (int j = 0; j < nc; j++) {
    for (int k = 0; k < n; k++) {
      for (int idx = all_p[k][j]; idx < all_p[k][j+1]; idx++) {
        int r = all_i[k][idx];
        if (mark[r] != j) {
          mark[r] = j;
          rows.push_back(r);
        }
      }
    }
    p_out[j+1] = p_out[j] + rows.size();
    // Clear marks for next column
    for (int r : rows) mark[r] = -1;
    rows.clear();
  }

  int nnz_out = p_out[nc];
  IntegerVector i_out(nnz_out);
  NumericVector  x_out(nnz_out);

  // ---- Pass 2: accumulate values ----
  int pos = 0;
  for (int j = 0; j < nc; j++) {
    for (int k = 0; k < n; k++) {
      for (int idx = all_p[k][j]; idx < all_p[k][j+1]; idx++) {
        int r = all_i[k][idx];
        if (mark[r] != j) {
          mark[r] = j;
          rows.push_back(r);
        }
        val[r] += all_x[k][idx];
      }
    }
    // dgCMatrix requires rows sorted within each column
    std::sort(rows.begin(), rows.end());
    for (int r : rows) {
      i_out[pos] = r;
      x_out[pos] = val[r];
      val[r] = 0.0;
      mark[r] = -1;
      pos++;
    }
    rows.clear();
  }

  // Build result dgCMatrix
  S4 res("dgCMatrix");
  res.slot("Dim") = IntegerVector::create(nr, nc);
  res.slot("i")   = i_out;
  res.slot("p")   = p_out;
  res.slot("x")   = x_out;

  // Preserve dimnames from the first matrix
  List dn0 = S4(matrices[0]).slot("Dimnames");
  res.slot("Dimnames") = dn0;

  return res;
}

// Per-LR communication layer kernel (computeCommunProb v2).
// Single fused pass over the CSC pattern of P.spatial:
//   v[k] = x0[k] * hill(L[row0[k]] * R[col0[k]]) * contact_mask[k]
//          [ * fAG[row]*fAG[col] * fAN[row]*fAN[col] ]   (agonist/antagonist factors)
// Multiplication grouping follows the baseline order exactly:
//   ((x0*h)*mask) * (AG_i*AG_j) * (AN_i*AN_j)  -> bitwise identical results.
// contact_mask: 1.0 at every position for diffusible LRs; adj.contact values
//   (1.0) at matched positions and 0.0 elsewhere for contact-dependent LRs
//   (masks the layer down to the contact pattern, as baseline `P1_Pspatial * adj.contact`).
// Returns list(i, p, x): compacted dgCMatrix slots (i 0-based; v != 0 kept; CSC order preserved).
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
// [[Rcpp::export]]
Rcpp::List cpp_prob_layer(Rcpp::NumericVector x0, Rcpp::IntegerVector row0,
                          Rcpp::IntegerVector col0,
                          Rcpp::NumericVector L, Rcpp::NumericVector R,
                          double Kh, double n,
                          Rcpp::NumericVector contact_mask,
                          Rcpp::Nullable<Rcpp::NumericVector> fAG,
                          Rcpp::Nullable<Rcpp::NumericVector> fAN,
                          int nC, int nthreads = 1) {
  R_xlen_t nnz = x0.length();
  if ((R_xlen_t)row0.length() != nnz || (R_xlen_t)col0.length() != nnz ||
      (R_xlen_t)contact_mask.length() != nnz)
    stop("x0, row0, col0, contact_mask must have the same length");
  bool ag = fAG.isNotNull(), an = fAN.isNotNull();
  Rcpp::NumericVector fag, fan;
  if (ag) fag = Rcpp::NumericVector(fAG.get());
  if (an) fan = Rcpp::NumericVector(fAN.get());

  Rcpp::NumericVector v(nnz);
  #pragma omp parallel for num_threads(nthreads) schedule(static)
  for (R_xlen_t k = 0; k < nnz; k++) {
    double lr = L[row0[k]] * R[col0[k]];
    double h;
    if (n == 1.0) h = lr / (Kh + lr);
    else { double ln = std::pow(lr, n), kn = std::pow(Kh, n); h = ln / (kn + ln); }
    v[k] = x0[k] * h * contact_mask[k];
  }
  if (ag || an) {
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (R_xlen_t k = 0; k < nnz; k++) {
      if (ag) v[k] = v[k] * (fag[row0[k]] * fag[col0[k]]);
      if (an) v[k] = v[k] * (fan[row0[k]] * fan[col0[k]]);
    }
  }

  // Compact (v != 0 kept); CSC order preserved, column counts -> p
  std::vector<int> cnt(nC, 0);
  for (R_xlen_t k = 0; k < nnz; k++)
    if (v[k] != 0.0) cnt[col0[k]]++;
  Rcpp::IntegerVector p(nC + 1);
  p[0] = 0;
  for (int j = 0; j < nC; j++) p[j + 1] = p[j] + cnt[j];
  Rcpp::IntegerVector i_out(p[nC]);
  Rcpp::NumericVector x_out(p[nC]);
  std::vector<int> off(p.begin(), p.end());
  for (R_xlen_t k = 0; k < nnz; k++) {
    if (v[k] != 0.0) {
      int pos = off[col0[k]]++;
      i_out[pos] = row0[k];
      x_out[pos] = v[k];
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = i_out,
                            Rcpp::Named("p") = p,
                            Rcpp::Named("x") = x_out);
}

// ============================================================================
// computeAvgCommunProb v2 (Plan A): gated group-level aggregation kernels.
//
// Replicates baseline computeAvgCommunProb_LR_{Avg,Sum} (modeling.R) element-for-element:
//   num[a,b] = sum_{i in a (sender row), j in b (receiver col)} v_ij
//              accumulated in CSC traversal order (j outer ascending, i inner ascending)
//   den[a,b] = count of STORED entries in the block (baseline binarized-crossprod;
//              explicit zeros count too, matching prob@x <- 1 + crossprod)
//   avg: val = num/den with 0/0 -> 0 (baseline NaN -> 0); sum: val = num
//   val *= indL[a]*indR[b]*indS[a]*indI[b]   (0/1 gate product = Prob_percent * cells.sr)
//
// Gate decisions are INTEGER threshold comparisons:
//   indL/indR: cnt >= thr[a] -- min.percent. thr[a] is precomputed in R against the
//     baseline's own mean()+signif() chain (decision is monotone non-decreasing in the
//     0/1 count cnt because group sizes are invariant under label permutation, so the
//     mean denominator n_a is a per-group constant; monotonicity is verified
//     exhaustively in tests_dev/test-computeAvgCommunProb.R).
//   indS/indI: sum >= min_cells_sr -- min.cells.sr. Sums of integer-valued doubles are
//     exact in double for any accumulation order.
// All gate quantities are integer-valued doubles: order cannot change them.
// The only order-sensitive quantity is num (~1 ulp vs the baseline's row-partial +
// BLAS dgemm chain; empirically verified in tests).
// Output layout: flat K*K with the SENDER group index fastest, f = a + b*K,
// i.e. exactly as.vector() of the K x K baseline matrix (column-major, a fastest).

#include <limits>
#include <algorithm>
#include <vector>

namespace {

void sc_check_group_args(const Rcpp::NumericVector& x, const Rcpp::IntegerVector& idx,
                         const Rcpp::IntegerVector& p, int nC, int K,
                         const Rcpp::NumericVector& suppL, const Rcpp::NumericVector& suppR,
                         const Rcpp::NumericVector& sr_out, const Rcpp::NumericVector& sr_in,
                         const Rcpp::IntegerVector& thr) {
  if ((R_xlen_t)idx.length() != x.length()) stop("i and x must have the same length");
  if (p.length() != (R_xlen_t)nC + 1) stop("p must have length nC + 1");
  if (suppL.length() != (R_xlen_t)nC || suppR.length() != (R_xlen_t)nC ||
      sr_out.length() != (R_xlen_t)nC || sr_in.length() != (R_xlen_t)nC)
    stop("suppL, suppR, sr_out, sr_in must have length nC");
  if (thr.length() != (R_xlen_t)K) stop("thr must have length k");
}

void sc_check_group_codes(const int* g, int nC, int K) {
  for (int c = 0; c < nC; c++)
    if (g[c] < 0 || g[c] >= K)
      stop("group codes must be 0-based integers in [0, k-1]");
}

// One gated aggregation pass for a single label assignment g (0-based codes).
void sc_boot_block_avg(const double* x, const int* idx, const int* p, int nC,
                       const double* suppL, const double* suppR,
                       const double* sr_out, const double* sr_in,
                       const int* g, int K,
                       const int* thr, double min_cells_sr, bool do_avg,
                       double* out) {
  std::vector<double> cntL(K, 0.0), cntR(K, 0.0), Sout(K, 0.0), Sin(K, 0.0);
  for (int c = 0; c < nC; c++) {
    const int a = g[c];
    cntL[a] += suppL[c];
    cntR[a] += suppR[c];
    Sout[a] += sr_out[c];
    Sin[a]  += sr_in[c];
  }
  std::vector<double> gL(K), gR(K), gS(K), gI(K);
  for (int a = 0; a < K; a++) {
    gL[a] = (cntL[a] >= (double)thr[a]) ? 1.0 : 0.0;
    gR[a] = (cntR[a] >= (double)thr[a]) ? 1.0 : 0.0;
    gS[a] = (Sout[a] >= min_cells_sr)   ? 1.0 : 0.0;
    gI[a] = (Sin[a]  >= min_cells_sr)   ? 1.0 : 0.0;
  }
  // Baseline early return when no group pair passes the min.percent gate
  // (sum(Prob_percent) == 0): result is all zeros regardless of the sr gate.
  bool any = false;
  for (int a = 0; a < K && !any; a++)
    for (int b = 0; b < K; b++)
      if (gL[a] * gR[b] * gS[a] * gI[b] != 0.0) { any = true; break; }
  if (!any) {
    std::fill(out, out + (std::size_t)K * K, 0.0);
    return;
  }
  std::vector<double> num((std::size_t)K * K, 0.0), den((std::size_t)K * K, 0.0);
  for (int j = 0; j < nC; j++) {
    const std::size_t gj = (std::size_t)g[j] * K;
    for (int t = p[j]; t < p[j + 1]; t++) {
      const std::size_t f = (std::size_t)g[idx[t]] + gj;
      num[f] += x[t];
      den[f] += 1.0;
    }
  }
  for (int b = 0; b < K; b++) {
    for (int a = 0; a < K; a++) {
      const std::size_t f = (std::size_t)a + (std::size_t)b * K;
      double val;
      if (do_avg) {
        val = den[f] > 0.0 ? num[f] / den[f]
                           : std::numeric_limits<double>::quiet_NaN();
        if (std::isnan(val)) val = 0.0;  // baseline: Prob.avg[is.nan(Prob.avg)] <- 0
      } else {
        val = num[f];
      }
      out[f] = val * (gL[a] * gR[b] * gS[a] * gI[b]);
    }
  }
}

}  // namespace

// Observed-labels aggregation for one LR layer (baseline averages section).
// group_int: 0-based group codes, length nC. avg_sum: 0 = "avg", 1 = "sum".
// [[Rcpp::export]]
Rcpp::NumericVector cpp_group_avg_obs(Rcpp::NumericVector x, Rcpp::IntegerVector idx,
                                      Rcpp::IntegerVector p,
                                      Rcpp::NumericVector suppL, Rcpp::NumericVector suppR,
                                      Rcpp::NumericVector sr_out, Rcpp::NumericVector sr_in,
                                      Rcpp::IntegerVector group_int, int K,
                                      Rcpp::IntegerVector thr, double min_cells_sr,
                                      int avg_sum, int nC) {
  sc_check_group_args(x, idx, p, nC, K, suppL, suppR, sr_out, sr_in, thr);
  if (group_int.length() != (R_xlen_t)nC) stop("group_int must have length nC");
  if (avg_sum != 0 && avg_sum != 1) stop("avg_sum must be 0 (avg) or 1 (sum)");
  sc_check_group_codes(group_int.begin(), nC, K);
  Rcpp::NumericVector out((R_xlen_t)K * K);
  sc_boot_block_avg(x.begin(), idx.begin(), p.begin(), nC,
                    suppL.begin(), suppR.begin(), sr_out.begin(), sr_in.begin(),
                    group_int.begin(), K, thr.begin(), min_cells_sr, avg_sum == 0,
                    out.begin());
  return out;
}

// Permutation aggregation for one LR layer (baseline permutation section).
// perm: nC x nboot integer matrix, column b = one sample.int(nC, nC) draw
// (1-based cell indices), i.e. the baseline `permutation` matrix passed as-is.
// Column b of the output is the k*k flat Pboot vector for boot b. Boots are
// independent (disjoint output columns, thread-local scratch), so results are
// bitwise identical for any nthreads.
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_group_avg_perm(Rcpp::NumericVector x, Rcpp::IntegerVector idx,
                                       Rcpp::IntegerVector p,
                                       Rcpp::NumericVector suppL, Rcpp::NumericVector suppR,
                                       Rcpp::NumericVector sr_out, Rcpp::NumericVector sr_in,
                                       Rcpp::IntegerVector group_int, int K,
                                       Rcpp::IntegerVector thr, double min_cells_sr,
                                       int avg_sum, Rcpp::IntegerMatrix perm,
                                       int nthreads = 1) {
  const int nC = perm.nrow();
  sc_check_group_args(x, idx, p, nC, K, suppL, suppR, sr_out, sr_in, thr);
  if (group_int.length() != (R_xlen_t)nC) stop("group_int must have length nC");
  if (avg_sum != 0 && avg_sum != 1) stop("avg_sum must be 0 (avg) or 1 (sum)");
  const int nboot = perm.ncol();
  Rcpp::NumericMatrix out((R_xlen_t)K * K, nboot);
  // perm 范围校验（并行区外）：条目必须落在 1..nC（sample.int(nC) 的值域）
  for (int c = 0; c < nC; c++)
    for (int b = 0; b < nboot; b++)
      if (perm(c, b) < 1 || perm(c, b) > nC)
        stop("perm entries must be cell indices in [1, nC]");
  const bool do_avg = avg_sum == 0;
  #pragma omp parallel for num_threads(nthreads) schedule(static)
  for (int b = 0; b < nboot; b++) {
    // 基线语义 gather：boot 标签 g_b[c] = group_int[perm(c, b)]（不是把 perm 当组码！）
    std::vector<int> g(nC);
    for (int c = 0; c < nC; c++) g[c] = group_int[perm(c, b) - 1];
    sc_boot_block_avg(x.begin(), idx.begin(), p.begin(), nC,
                      suppL.begin(), suppR.begin(), sr_out.begin(), sr_in.begin(),
                      g.data(), K, thr.begin(), min_cells_sr, do_avg,
                      &out(0, b));
  }
  return out;
}
