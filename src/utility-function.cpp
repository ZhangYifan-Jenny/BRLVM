#include "utility-function.h"

// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd rwish1(double v, Eigen::MatrixXd S) {
  arma::mat armaS = arma::mat(S.data(), S.rows(), S.cols(), false, true);
  
  int p = armaS.n_rows;
  arma::mat CC = arma::chol(armaS,"lower");  // Cholesky decomposition
  arma::mat Z = arma::zeros(p, p);   // Initialize Z matrix
  
  // Set the diagonal elements from a chi-squared distribution
  for (int i = 0; i < p; ++i) {
    Z(i, i) = std::sqrt(arma::chi2rnd(v - i)); // Chi-squared distribution
  }
  
  // Fill upper triangular elements with normal random variables
  //for (int i = 0; i < p; ++i) {
  //  for (int j = i + 1; j < p; ++j) {
  //    Z(i, j) = arma::randn();
  //  }
  //}
  arma::vec sampling = arma::randn((p * (p - 1)) / 2);
  
  int k = 0; // Index for sampling vector
  
  // Assign the sampling to Z lower matrix using a loop
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      Z(i, j) = sampling(k);
      ++k; // Increment the index for sampling vector
    }
  }
  arma::mat W = CC * Z * Z.t() * CC.t();
  Eigen::MatrixXd result = Eigen::Map<Eigen::MatrixXd>(W.memptr(), W.n_rows, W.n_cols);
 
  return  result;
}
Rcpp::NumericVector rinvgauss(int n, Rcpp::NumericVector mu, double lambda) {
  //arma::vec mu_vec(mu.data(), mu.size(), false);
  arma::vec mu_vec(mu.begin(), mu.size(), false);
  arma::vec un = arma::randu<arma::vec>(n);
  arma::vec Xi = arma::chi2rnd<arma::vec>(arma::ones<arma::vec>(n) * 1);
  arma::vec f = mu_vec / (2 * lambda) % (2 * lambda * arma::ones<arma::vec>(mu_vec.size()) + mu_vec % Xi + arma::sqrt(4 * lambda * mu_vec % Xi + arma::square(mu_vec) % arma::square(Xi)));
  arma::vec s = pow(mu_vec, 2) / f;
  arma::vec result_vec = arma::zeros<arma::vec>(n);
  arma::uvec indices = arma::find(un < mu_vec / (mu_vec + s));
  result_vec(indices) = s(indices);
  result_vec(arma::find(un >= mu_vec / (mu_vec + s))) = f(arma::find(un >= mu_vec / (mu_vec + s)));
  //Rcpp::NumericVector result(result_vec.begin(), result_vec.end());
  Rcpp::NumericVector result = Rcpp::wrap(result_vec);
  return result;
}

//X[:,sub2[i,:]]
Eigen::MatrixXd  selectColumns (Eigen::MatrixXd X, Eigen::MatrixXi sub1, int i) {
  // Count the number of True values in the i-th row of sub1
  int numTrue = (sub1.row(i).array() == 1).count() ;
  // Create a new matrix to store selected columns
  Eigen::MatrixXd selected(X.rows(), numTrue);                                                                                   
  // Iterate over columns and copy the ones where sub1[i, col] is true
  int colIndex = 0;
  for (int col = 0; col < X.cols(); ++col) {
    if (sub1(i, col) == 1) {
      selected.col(colIndex++) = X.col(col);
    }
  }
  return selected;
};

//C[rf, :][:,rf] : submatrix(C,rf,rf)
//choosing row and col from matrix C
//want all rows: row=Eigen::VectorXi::Ones(1)
Eigen::MatrixXd  submatrix (Eigen::MatrixXd C, Eigen::VectorXi row, Eigen::VectorXi col) {
  if(row.size() != C.rows()) {
    row = Eigen::VectorXi::Constant(C.rows(), row[0]); 
  }
  if(col.size() != C.cols()) {
    col = Eigen::VectorXi::Constant(C.cols(), col[0]); 
  }
  Eigen::MatrixXd subC(row.sum(),col.sum());
  int rowIndex = 0;
  for(int i = 0; i < row.sum(); ++i) {
    if(row[i]) {
      int colIndex = 0;
      for(int j = 0; j < col.sum(); ++j) {
        if(col[j]) {
          subC(i,j) = C(rowIndex, colIndex++);
        }
      }
      rowIndex++;
    }
  }
  
  return subC;
};

Eigen::MatrixXd hpd_interval(Eigen::MatrixXd mc_samples, double prob ) {
  int n = mc_samples.rows();
  int p = mc_samples.cols();
  int m = int(prob * n);
  Eigen::MatrixXd intervals(p, 2);
  
  for (int i = 0; i < p; ++i) {
    Eigen::VectorXd samples_col = mc_samples.col(i);
    std::sort(samples_col.data(), samples_col.data() + samples_col.size());
    
    double min_width = std::numeric_limits<double>::max();
    int hpd_index = 0;
    
    for (int j = 0; j < n - m; ++j) {
      double width = samples_col(j + m) - samples_col(j);
      if (width < min_width) {
        min_width = width;
        hpd_index = j;
      }
    }
    
    intervals(i, 0) = samples_col(hpd_index);
    intervals(i, 1) = samples_col(hpd_index + m);
  }
  return intervals;
}

Eigen::VectorXd compute_psrf(Eigen::MatrixXd chain) { //chain is passed as an Eigen::MatrixXd
  int N = chain.rows();
  int num_vars = chain.cols();
  Eigen::VectorXd r_hats(num_vars);
  
  for (int i = 0; i < num_vars; ++i) {
    Eigen::VectorXd variable_chain = chain.col(i);
    
    // Split the chain into two halves
    Eigen::VectorXd chain1 = variable_chain.head(N / 2);
    Eigen::VectorXd chain2 = variable_chain.tail(N / 2);
    
    // Calculate the within-chain variance W
    double W = 0.5 * ((chain1.array() - chain1.mean()).square().sum() / (chain1.size() - 1) +
                      (chain2.array() - chain2.mean()).square().sum() / (chain2.size() - 1));
    
    double r_hat = 0;
    if (W != 0) {
      // Calculate the between-chain variance B
      double mean_diff = chain1.mean() - chain2.mean();
      double B = (N / 2.0) * std::pow(mean_diff, 2);
      // Compute the total variance V
      double V = (1 - 1.0 / N) * W + (1.0 / N) * B;
      // Calculate the R-hat (PSRF) for the current variable
      r_hat = std::sqrt(V / W);
    }
    
    r_hats[i] = r_hat;
  }
  
  return r_hats;
}