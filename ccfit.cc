#include "ccfit.hh"

#include <iostream>

#include <Eigen/Core>

using namespace Eigen;

ClosedCurve closedCurveFit(const std::vector<double> &points,
                           const std::vector<size_t> &important_indices,
                           size_t degree, size_t n_segments) {
  size_t n = points.size() / 3;
  Map<const Matrix<double,Dynamic,3,RowMajor>> b(&points[0], n, 3);

  ClosedCurve result;
  result.degree = degree;

  // TODO

  return result;
}

