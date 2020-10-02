#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ccfit.hh"

using PointSet = std::vector<double>;

std::vector<PointSet> readParams(std::string filename) {
  const std::array<std::string, 3> types = { "top", "middle", "bottom" };
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<PointSet> result(types.size());
  std::string line, word;
  while (!f.eof()) {
    std::getline(f, line);
    f >> std::ws;
    if (line.empty())
      continue;
    std::istringstream ss(line);
    ss >> word;
    auto it = std::find(types.begin(), types.end(), word);
    if (it == types.end())
      continue;
    size_t index = std::distance(types.begin(), it);
    std::array<double, 3> p;
    ss >> p[0] >> p[1] >> p[2];
    result[index].insert(result[index].end(), p.begin(), p.end());
  }
  return result;
}

void writeCurve(const ClosedCurve &curve, std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  // Sampling
  double range = curve.knots.back() - curve.knots.front();
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / resolution; // does not reach 1
    u = curve.knots.front() + range * u;
    auto p = curve.eval(u);
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f << "l ";
  for (size_t i = 1; i <= resolution; ++i)
    f << i << ' ';
  f << 1 << std::endl;

  // Control polygon
  // for (size_t i = 0; i < curve.cpts.size(); i += 3)
  //   f << "v " << curve.cpts[i] << ' ' << curve.cpts[i+1] << ' ' << curve.cpts[i+2] << std::endl;
  // f << "l ";
  // for (size_t i = 1; i <= curve.cpts.size() / 3; ++i)
  //   f << resolution + i << ' ';
  // f << resolution + 1 << std::endl;
}

void writePoints(const std::vector<PointSet> &points, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t count = 0;
  for (const auto &ps : points)
    for (size_t i = 0; i < ps.size(); i += 3) {
      f << "v " << ps[i] << ' ' << ps[i+1] << ' ' << ps[i+2] << std::endl;
      ++count;
    }
  for (size_t i = 1; i <= count; ++i)
    f << "p " << i << std::endl;
}

void writeSurface(const std::vector<ClosedCurve> &curves, std::string filename, size_t resolution) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &curve : curves) {
    double range = curve.knots.back() - curve.knots.front();
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / resolution; // does not reach 1
      u = curve.knots.front() + range * u;
      auto p = curve.eval(u);
      f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }
  }

  // Top
  f << "f";
  for (size_t i = 1; i <= resolution; ++i)
    f << ' ' << i;
  f << ' ' << 1 << std::endl;

  // Bottom
  f << "f";
  for (size_t i = resolution * (curves.size() - 1) + 1; i <= resolution * curves.size(); ++i)
    f << ' ' << i;
  f << ' ' << resolution * (curves.size() - 1) + 1 << std::endl;

  // Side
  for (size_t level = 1; level < curves.size(); ++level)
    for (size_t i = 1; i <= resolution; ++i) {
      size_t ip = (i % resolution) + 1;
      size_t j = (level - 1) * resolution + i;
      size_t jp = (level - 1) * resolution + ip;
      size_t k = level * resolution + i;
      size_t kp = level * resolution + ip;
      f << "f " << j << ' ' << jp << ' ' << kp << ' ' << k << std::endl;
    }
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.params> [# of segments]" << std::endl;
    return 1;
  }

  size_t n_segments = 8;
  if (argc == 3)
    n_segments = std::atoi(argv[2]);

  auto params = readParams(argv[1]);
  writePoints(params, "/tmp/points.obj");

  std::vector<ClosedCurve> curves;
  for (size_t i = 0; i < 3; ++i) {
    const auto &points = params[i];
    auto curve = closedCurveFit(points, {8, 17}, 3, n_segments);
    curves.push_back(curve);
    writeCurve(curve, std::string("/tmp/curve-") + std::to_string(i) + ".obj", 100);
  }

  writeSurface(curves, "/tmp/surface.obj", 100);
}
