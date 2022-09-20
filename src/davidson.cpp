#include "src/davidson.hpp"

namespace davidson {

template <>
class DavidsonEigenSolver<Eigen::MatrixXd>;

template <>
class DavidsonEigenSolver<Eigen::MatrixXcd>;

}  // namespace davidson