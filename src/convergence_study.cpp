#define STAMP 0
#include "integrators.h"
#include "vectors.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

Vec3 E(Vec3 x, double t);
Vec3 B(Vec3 x, double t);

int main(int argc, char *argv[]) {
  std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(7)
            << std::right;
  std::ofstream solution;
  solution.open("../output/one_conv.dat");

  double q[] = {1.};
  double m[] = {1.};
  Vec3 xrk[] = {Vec3{0.3, 0., 0.}};
  Vec3 vrk[] = {Vec3{0., 1., 0.}};
  Vec3 xbo[] = {Vec3{0.3, 0., 0.}};
  Vec3 vbo[] = {Vec3{0., 1., 0.}};
  int numpoints = 10000;
  double dt = 0.03;
  double t = 0.;
#if STAMP
  std::cout << t << " " << std::setw(15) << xrk[0].x << " " << std::setw(15)
            << xrk[0].y << " " << std::setw(15) << xrk[0].z << " "
            << std::setw(15) << xbo[0].x << " " << std::setw(15) << xbo[0].y
            << " " << std::setw(15) << xrk[0].z << std::endl;
#endif
  solution << t << " " << std::setw(15) << xrk[0].x << " " << std::setw(15)
           << xrk[0].y << " " << std::setw(15) << xrk[0].z << " "
           << std::setw(15) << xbo[0].x << " " << std::setw(15) << xbo[0].y
           << " " << std::setw(15) << xrk[0].z << "\n";

  for (int i = 0; i < numpoints; i++) {
    RK4_step(t, xrk, vrk, q, m, E, B, dt, 1, xrk, vrk);
    boris_step(t, xbo, vbo, q, m, E, B, dt, 1, xbo, vbo);

    t += dt;
#if STAMP
    std::cout << t << " " << std::setw(15) << xrk[0].x << " " << std::setw(15)
              << xrk[0].y << " " << std::setw(15) << xrk[0].z << " "
              << std::setw(15) << xbo[0].x << " " << std::setw(15) << xbo[0].y
              << " " << std::setw(15) << xrk[0].z << std::endl;
#endif
    solution << t << " " << std::setw(15) << xrk[0].x << " " << std::setw(15)
             << xrk[0].y << " " << std::setw(15) << xrk[0].z << " "
             << std::setw(15) << xbo[0].x << " " << std::setw(15) << xbo[0].y
             << " " << std::setw(15) << xrk[0].z << "\n";
  }

  solution.close();

  return 0;
}

// Vec3 E(Vec3 x, double t) { return Vec3{0., 0.5, 0.}; }
//
// Vec3 B(Vec3 x, double t) { return Vec3{0., 0., 1.}; }

Vec3 B(Vec3 x, double t) { return Vec3{0., 0., sqrt(x.x * x.x + x.y * x.y)};
}

Vec3 E(Vec3 x, double t) {
  double d = sqrt(x.x * x.x + x.y * x.y);
  d = 1.e-2 / (d * d * d);
  return Vec3{-x.x * d, -x.y * d, 0.};
}
