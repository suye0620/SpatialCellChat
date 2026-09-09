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
