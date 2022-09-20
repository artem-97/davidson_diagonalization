#pragma once

#include <gtest/gtest.h>

#include <iostream>

#include "src/davidson.hpp"
#include "tests/asserts.hpp"

class DavidsonTest : public ::testing::Test {};

TEST_F(DavidsonTest, davdiag_real) {
  constexpr auto n{1200};
  constexpr auto tol{1e-8};
  constexpr auto sparsity{0.0001};

  constexpr auto k{8};
  constexpr auto eig{4};
  Eigen::MatrixXd A(n, n);
  A.setZero();
  for (size_t i = 0; i < n; ++i) A(i, i) = i + 1;

  A += sparsity * Eigen::MatrixXd::Random(n, n);
  A = 0.5 * (A.transpose() + A);

  davidson::DavidsonEigenSolver dav_es(A, k, eig);
  Eigen::SelfAdjointEigenSolver<decltype(A)> eig_es(A);
  for (size_t i = 0; i < eig; ++i)
    ASSERT_APPROX_EQ(dav_es.eigenvalues()[i], eig_es.eigenvalues()[i], tol);
}

TEST_F(DavidsonTest, davdiag_cmplx) {
  constexpr auto n{200};
  constexpr auto tol{1e-8};
  constexpr auto sparsity{0.0001};

  constexpr auto k{8};
  constexpr auto eig{4};
  Eigen::MatrixXcd A(n, n);
  A.setZero();
  for (size_t i = 0; i < n; ++i) A(i, i) = i + 1;

  A += sparsity * Eigen::MatrixXcd::Random(n, n);
  A = 0.5 * (A.transpose() + A);

  davidson::DavidsonEigenSolver dav_es(A, k, eig);
  Eigen::SelfAdjointEigenSolver<decltype(A)> eig_es(A);
  for (size_t i = 0; i < eig; ++i) {
    auto x = dav_es.eigenvalues()[i];
    auto y = dav_es.eigenvalues()[i];
    ASSERT_APPROX_EQ(std::abs(x), std::abs(y), tol);
  }
}