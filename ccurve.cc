#include "ccurve.hh"

#include <Eigen/Core>

using namespace Eigen;

std::array<double, 3> ClosedCurve::eval(double) const {
  Map<const Matrix<double,Dynamic,3,RowMajor>> cp(&cpts[0], cpts.size() / 3, 3);
  Vector3d p;

  // TODO

  std::array<double, 3> result;
  std::copy_n(p.data(), 3, result.begin());
  return result;
}

