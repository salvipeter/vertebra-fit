#pragma once

#include <array>
#include <vector>

// `knots` has n + 1 values
// `cpts` has (3*)n values - it does not include the repeated first control point at the end
struct ClosedCurve {
  size_t degree;
  std::vector<double> knots;
  std::vector<double> cpts; // x1 y1 z1 x2 y2 z2 ...

  size_t span(double u) const;
  std::vector<double> basis(size_t span, double u) const;
  std::array<double, 3> eval(double u) const;
};
