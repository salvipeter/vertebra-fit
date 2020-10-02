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
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / resolution; // does not reach 1
    auto p = curve.eval(u);
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f << "l ";
  for (size_t i = 1; i <= resolution; ++i)
    f << i << ' ';
  f << '1' << std::endl;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <model.params>" << std::endl;
    return 1;
  }

  auto params = readParams(argv[1]);
  for (size_t i = 0; i < 3; ++i) {
    const auto &points = params[i];
    auto curve = closedCurveFit(points, {0, 1});
    writeCurve(curve, std::string("/tmp/") + std::to_string(i) + ".obj", 100);
  }
}
