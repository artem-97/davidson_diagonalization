#include "src/davidson.hpp"
#include "tests/davidson.hpp"
#include "tests/qr.hpp"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// https://joshuagoings.com/2013/08/23/davidsons-method/
/*
int main() {
  //   constexpr auto n{1200};
  constexpr auto n{6};

  constexpr auto tol{1e-8};
  constexpr auto mmax{n / 2};
  constexpr auto sparsity{0.0001};
  Eigen::MatrixXd A(n, n);
  A.setZero();
  for (size_t i = 0; i < n; ++i) A(i, i) = i + 1;
  A(2, 3) = 0.004;
  A(3, 2) = 0.004;

  A += sparsity * Eigen::MatrixXd::Random(n, n);
  A = 0.5 * (A.transpose() + A);


  //   auto k = 2000;
  //   A = Eigen::MatrixXd::Random(k, k);
  //   A = 0.5 * (A + A.transpose());
  //   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  //   auto e = es.eigenvalues();

  davidson::DavidsonEigenSolver es(A, 2, 1);
  //   auto [q, r] = davidson::DavidsonEigenSolver::QR(A);
  //   std::cout << q << std::endl;
  //   std::cout << "----\n";
  //   std::cout << r << std::endl;

  //   std::cout << A.col(0).normalized() << std::endl;
  //   std::cout << davidson::DavidsonEigenSolver::eye(4, 8) << std::endl;
}
*/