#include "integrators.h"
#include "vectors.h"

#include <cmath>

Vec3 E(Vec3 x, double t);
Vec3 B(Vec3 x, double t);

int main(int argc, char *argv[]) { 
  double dt = 1.;
  int n_iter = 500;

  return 0;
}

Vec3 E(Vec3 x, double t) { return Vec3{0., 0., sqrt(x.x * x.x + x.y * x.y)}; }

Vec3 B(Vec3 x, double t) {
  double d = sqrt(x.x * x.x + x.y * x.y);
  d = 1. / (d * d * d);
  return Vec3{-x.x * d, -x.y * d, 0.};
}
