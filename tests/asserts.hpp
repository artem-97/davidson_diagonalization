#pragma once

#include <gtest/gtest.h>

#define ASSERT_APPROX_EQ(x, y, tol) ASSERT_LT(std::abs(x - y), tol)

#define ASSERT_MATRIX_EQ(M1, M2)                   \
  ASSERT_EQ(M1.rows(), M2.rows());                 \
  ASSERT_EQ(M1.cols(), M2.cols());                 \
  for (int i = 0; i < M1.rows(); ++i) {            \
    for (int j = 0; j < M1.cols(); ++j) {          \
      ASSERT_APPROX_EQ(M1(i, j), M2(i, j), 1e-10); \
    }                                              \
  }
