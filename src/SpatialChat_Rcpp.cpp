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
