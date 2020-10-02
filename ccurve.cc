#include "ccurve.hh"

#include <Eigen/Core>

using namespace Eigen;

size_t ClosedCurve::span(double u) const {
  if (u <= knots.front())
    return 0;
  if (u >= knots.back())
    return knots.size() - 1;
  auto it = std::upper_bound(knots.begin(), knots.end(), u);
  return std::distance(knots.begin(), it) - 1;
}

std::vector<double> ClosedCurve::basis(size_t span, double u) const {
  size_t n = knots.size() - 1;
  double range = knots.back() - knots.front();
  std::vector<double> coeff, left(degree + 1), right(degree + 1);
  coeff.push_back(1);
  for(size_t j = 1; j <= degree; ++j) {
    if (span + 1 < j)
      left[j] = range + (u - knots.front()) - knots[span+1+n-j];
    else
      left[j] = u - knots[span+1-j];
    if (span + j >= n)
      right[j] = range + (knots[span+j-n] - knots.front()) - u;
    else
      right[j] = knots[span+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
  return coeff;
}

std::array<double, 3> ClosedCurve::eval(double u) const {
  size_t n_cp = cpts.size() / 3;
  Map<const Matrix<double,Dynamic,3,RowMajor>> cp(&cpts[0], n_cp, 3);

  size_t s = span(u);
  auto coeff = basis(s, u);
  Vector3d p(0, 0, 0);
  for(size_t i = 0; i <= degree; ++i)
    if (s + i < degree)
      p += cp.row(n_cp + s - degree + i) * coeff[i];
    else
      p += cp.row((s - degree + i) % n_cp) * coeff[i];

  std::array<double, 3> result;
  std::copy_n(p.data(), 3, result.begin());
  return result;
}

