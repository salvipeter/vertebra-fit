#pragma once

#include "ccurve.hh"

ClosedCurve closedCurveFit(const std::vector<double> &points, // x1 y1 z1 x2 y2 z2 ...
                           const std::vector<size_t> &important_indices,
                           size_t degree, size_t n_cp);
