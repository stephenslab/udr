// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

const double LOG_2PI = std::log(2.0 * M_PI);

// Softmax functions: yi = exp(xi) / sum(exp(xj))
inline arma::vec
softmax(const arma::vec & x)
{
    // Calculate exp()
    // Subtract the max - this prevents overflow, which happens for x ~ 1000
    arma::vec y = arma::exp(x - arma::max(x));
    // Renormalise
    y /= arma::sum(y);
    return y;
}

// function for "shrinking" the covariance matrix, to get $\hat U_k$.
inline arma::mat
shrink_cov(const arma::mat & V, const double & eps)
{
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, V);
    for (arma::uword i = 0; i < eigval.n_elem; ++i) {
        eigval(i) = (eigval(i) > 1.0) ? eigval(i) : (1.0 + eps);
    }
    return eigvec * diagmat(eigval) * trans(eigvec);
}

inline arma::vec
dmvnorm_mat(const arma::mat & x,
  const arma::vec           & mean,
  const arma::mat           & sigma,
  bool                      logd     = false,
  bool                      inversed = false)
{
    double xdim = static_cast<double>(x.n_rows);

    arma::vec out(x.n_cols);
    arma::mat rooti;

    // we have previously computed rooti
    // in R eg rooti <- backsolve(chol(sigma), diag(ncol(x)))
    if (inversed) { rooti = sigma; } else {
        try {
            rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
        } catch (const std::runtime_error & error) {
            if (logd) out.fill(-arma::datum::inf);
            else out.fill(0.0);
            for (arma::uword i = 0; i < x.n_cols; ++i)
                if (arma::accu(arma::abs(x.col(i) - mean)) < 1e-6) out.at(i) = arma::datum::inf;
            return out;
        }
    }
    double rootisum  = arma::sum(arma::log(rooti.diag()));
    double constants = -(xdim / 2.0) * LOG_2PI;

    for (unsigned i = 0; i < x.n_cols; i++) {
        arma::vec z = rooti * (x.col(i) - mean);
        out.at(i) = constants - 0.5 * arma::sum(z % z) + rootisum;
    }

    if (logd == false) {
        out = arma::exp(out);
    }
    return out;
}

// Truncated Eigenvalue Extreme deconvolution
class TEEM
{
public:
    TEEM(const arma::mat & X_mat,
      const arma::vec    & w_vec,
      const arma::cube   & U_cube) :
        X_mat(X_mat), w_vec(w_vec)
    {
        T_cube = U_cube;
        for (unsigned j = 0; j < T_cube.n_slices; ++j) {
            T_cube.slice(j) += arma::eye(arma::size(T_cube.slice(j)));
        }
    }

    ~TEEM(){ }

    arma::vec
    get_objective(){ return objective; }

    arma::vec
    get_maxd(){ return maxd; }

    arma::vec
    get_w(){ return w_vec; }

    arma::cube
    get_U()
    {
        arma::cube U_cube = T_cube;
        for (unsigned j = 0; j < U_cube.n_slices; ++j) {
            U_cube.slice(j) -= arma::eye(arma::size(U_cube.slice(j)));
        }
        return U_cube;
    }

    int
    fit(const int & maxiter, const double & converge_tol, const double & eigen_tol, const bool & verbose)
    {
        // initialize to store progress
        objective.zeros(maxiter);
        maxd.zeros(maxiter);
        int iter_out;

        // Get the number of samples (n) and the number of mixture components (k)
        int n = X_mat.n_rows;
        // int m = X_mat.n_cols;
        int k = w_vec.size();

        for (unsigned iter = 0; iter < maxiter; ++iter) {
            // store parameters and likelihood in the previous step
            arma::vec w0_vec = w_vec;

            // E-step: calculate posterior probabilities using the current mu and sigmas
            // arma::mat logP = mat(n, k, arma::fill::zeros); // n by k matrix;
            arma::mat logP = arma::zeros<arma::mat>(n, k); // n by k matrix
            for (unsigned j = 0; j < k; ++j) {
                logP.col(j) = log(w_vec(j)) + dmvnorm_mat(trans(X_mat), arma::zeros<arma::vec>(
                        X_mat.n_cols), T_cube.slice(j), true); // ??
            }
            // softmax for renormalization
            // arma::mat P_mat = mat(k, n, arma::fill::zeros); // k by n matrix. because of row/col vec converting
            arma::mat P_mat = arma::zeros<arma::mat>(k, n); // k by n matrix. because of row/col vec converting

            for (arma::uword i = 0; i < n; ++i) {
                arma::colvec y = arma::conv_to<arma::colvec>::from(logP.row(i));
                P_mat.col(i) = softmax(y);
            }
            P_mat = trans(P_mat); // n by k matrix

            // M-step:
            for (unsigned j = 0; j < k; ++j) {
                T_cube.slice(j) = trans(X_mat) * (P_mat.col(j) % X_mat.each_col()) / accu(P_mat.col(j));
                T_cube.slice(j) = shrink_cov(T_cube.slice(j), eigen_tol);
            }
            // update mixture weights
            w_vec = arma::conv_to<arma::colvec>::from(sum(P_mat, 0)) / n; // 0:sum by column;

            // Compute log-likelihood at the current estimates
            double f = compute_loglik();

            // Check stopping criterion
            double d = max(abs(w_vec - w0_vec));
            maxd(iter)      = d;
            objective(iter) = f;
            iter_out        = iter;

            if (d < converge_tol) {
                break;
            }
        }
        objective.resize(iter_out + 1);
        maxd.resize(iter_out + 1);
        return 0;
    } // fit

private:
    arma::mat X_mat;
    arma::vec w_vec;
    arma::cube T_cube;
    arma::vec objective;
    arma::vec maxd;
    double
    compute_loglik()
    {
        int n = X_mat.n_rows;
        int k = w_vec.size();

        arma::vec y = arma::zeros<arma::vec>(n);
        for (unsigned j = 0; j < k; ++j) {
            y = y + w_vec(j) * dmvnorm_mat(trans(X_mat), arma::zeros<arma::vec>(X_mat.n_cols), T_cube.slice(j));
        }
        return (sum(log(y)));
    }
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List
fit_teem_rcpp(Rcpp::NumericMatrix X_mat,
  Rcpp::NumericVector             w_vec,
  Rcpp::NumericVector             U_3d,
  int                             maxiter,
  double                          converge_tol,
  double                          eigen_tol,
  bool                            verbose)
{
    // Convert R 3d array to Rcpp cube
    if (Rf_isNull(U_3d.attr("dim"))) {
        throw std::invalid_argument(
                  "U_3d has to be a 3D array");
    }
    Rcpp::IntegerVector dimU = U_3d.attr("dim");
    arma::cube U_cube(U_3d.begin(), dimU[0], dimU[1], dimU[2]);
    //
    TEEM teem(Rcpp::as<arma::mat>(X_mat),
      Rcpp::as<arma::vec>(w_vec),
      U_cube);
    teem.fit(maxiter, converge_tol, eigen_tol, verbose);
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named("w")         = teem.get_w(),
        Rcpp::Named("U")         = teem.get_U(),
        Rcpp::Named("objective") = teem.get_objective(),
        Rcpp::Named("maxd")      = teem.get_maxd());
    return res;
}
