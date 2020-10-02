#pragma once

#include <array>
#include <vector>

struct ClosedCurve {
  size_t degree;
  std::vector<double> knots;
  std::vector<double> cpts; // x1 y1 z1 x2 y2 z2 ...

  std::array<double, 3> eval(double) const;
};
