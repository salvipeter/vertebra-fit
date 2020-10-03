#include "ccfit.hh"

#include <cassert>

#include <Eigen/Core>
#include <Eigen/LU>

using namespace Eigen;

ClosedCurve closedCurveFit(const std::vector<double> &points,
                           const std::vector<size_t> &important_indices,
                           size_t degree, size_t n_cp) {
  size_t m = points.size() / 3;
  assert(n_cp <= m);
  Map<const Matrix<double,Dynamic,3,RowMajor>> data(&points[0], m, 3);
  MatrixXd b = data;

  ClosedCurve result;
  result.degree = degree;

  std::vector<double> u;
  u.push_back(0);
  double total = 0;
  for (size_t i = 1; i < m; ++i) {
    u.push_back(std::sqrt((b.row(i) - b.row(i - 1)).norm()));
    total += u.back();
  }
  u.push_back(std::sqrt((b.row(0) - b.row(m - 1)).norm()));
  total += u.back();
  for (size_t i = 1; i < m; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1;

  result.knots.push_back(0);
  for (size_t i = 1; i < n_cp; ++i) {
    double d = (double)m / n_cp; // >= 1
    size_t j = d * i;            // >= 1
    double alpha = d * i - j;
    double knot = (1 - alpha) * u[j-1] + alpha * u[j];
    result.knots.push_back(knot);
  }
  result.knots.push_back(1);

  u.pop_back();

  MatrixXd A = MatrixXd::Zero(m, n_cp);
  for (size_t i = 0; i < m; ++i) {
    size_t s = result.span(u[i]);
    auto coeff = result.basis(s, u[i]);
    for (size_t j = 0; j <= degree; ++j)
      if (s + j < degree)
        A(i, n_cp + s - degree + j) = coeff[j];
      else
        A(i, (s - degree + j) % n_cp) = coeff[j];
  }

  double large_number = 1000;
  for (size_t i : important_indices) {
    A.row(i) *= large_number;
    b.row(i) *= large_number;
  }

  MatrixXd x = A.fullPivLu().solve(b);
  for (size_t i = 0; i < n_cp; ++i)
    for (size_t j = 0; j < 3; ++j)
      result.cpts.push_back(x(i, j));

  return result;
}

