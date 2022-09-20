#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace davidson {

namespace detail {

template <typename MatrixT>
std::pair<MatrixT, MatrixT> QR(MatrixT A) {
  Eigen::HouseholderQR<MatrixT> qr(A);
  auto m = A.rows();
  auto n = A.cols();
  MatrixT Q = qr.householderQ();
  if (m > n) {
    Q = Q * MatrixT::Identity(m, n);
  }
  MatrixT R = Q.transpose() * A;
  return {Q, R};
}

}  // namespace detail

template <typename MatrixT>
class DavidsonEigenSolver {
 public:
  DavidsonEigenSolver(MatrixT M, size_t k, size_t eig)
      : M_(std::move(M)),
        n_(M_.rows()),
        mmax_(M_.rows() / 2),
        k_(k),
        eig_(eig) {
    assert(k_ >= eig_);
    run();
  }
  inline auto eigenvalues() { return eigenvalues_; }

 private:
  static MatrixT zeros(size_t n, size_t m) {
    MatrixT E(n, m);
    E.setZero();
    return E;
  }

  static MatrixT eye(size_t n, size_t m) {
    auto E = zeros(n, m);
    for (size_t i = 0; i < std::min(n, m); ++i) E(i, i) = 1;
    return E;
  }

  inline static MatrixT eye(size_t n) { return eye(n, n); }

  void run() {
    auto t = eye(n_, k_);
    auto V = zeros(n_, n_);
    auto I = eye(n_);

    Eigen::VectorXd theta_old, theta;

    for (size_t m = k_; m < mmax_; m += k_) {
      if (m <= k_) {
        for (size_t j = 0; j < k_; ++j) V.col(j) = t.col(j).normalized();
        theta_old = Eigen::VectorXd::Ones(eig_, 1);
      } else if (m > k_)
        theta_old = theta.head(eig_);

      MatrixT Vm = V.leftCols(m);
      auto [Q, R] = detail::QR(Vm);
      V.leftCols(m) = Q;
      auto T = Q.transpose() * M_ * Q;

      Eigen::SelfAdjointEigenSolver<MatrixT> es(T);
      theta = es.eigenvalues();
      auto s = es.eigenvectors();

      // std::cout << "---------------" << std::endl;
      // std::cout << theta[0] << " " << theta[1] << std::endl;
      // std::cout << "Vcols" << V.cols() << " :" << m << std::endl;
      for (size_t j = 0; j < k_; ++j) {
        auto fact = V.leftCols(m) * s.col(j);
        auto w = (M_ - theta[j] * I) * fact;
        // auto d = theta[j] - M_(j, j);
        // std::cout << "W" << w.cols() << " " << w.rows() << std::endl;
        auto q = w / (theta(j) - M_(j, j));
        // std::cout << "Q" << q.rows() << " " << q.cols() << std::endl;

        V.col(m + j) = q;
      }
      // break;
      auto norm = (theta.head(eig_ - 1) - theta_old.head(eig_ - 1)).norm();
      if (norm < tol_) break;
    }

    eigenvalues_.resize(theta.size());
    for (int i = 0; i < theta.size(); ++i) {
      eigenvalues_[i] = theta[i];
    }
  }

 private:
  MatrixT M_;
  size_t n_;
  size_t mmax_;
  double tol_{1e-8};
  size_t k_;  // initial guess vectors
  size_t eig_;

  std::vector<double> eigenvalues_;
};

}  // namespace davidson