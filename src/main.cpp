#include "conversions.h"
#include "init.h"
#include "integrators.h"
#include "utils.h"
#include "vectors.h"

#include <fstream>
#include <iomanip>
#include <iostream>

// TODO: Set up system for time limit
int main() {
  constexpr double dt = time_to_sim_units(STEP_SIZE);
  double *q = new double[N_PART];
  double *m = new double[N_PART];
  Vec3 *x = new Vec3[N_PART];
  Vec3 *v = new Vec3[N_PART];
  initialize_starting_state(x, v, q, m);

  std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(7)
            << std::right;
  std::ofstream solution;
  // TODO: Organize file output
  solution.open("../output/main.dat");

  log_units();

  double t = 0.;
  log_advancement(-1, STEPS_LIMIT);

  // TODO: Fix for N particles
  Vec3 x_cgs = space_to_cgs_units(x[0]);
  solution << std::setw(6) << time_to_cgs_units(t) << " " << std::setw(15)
           << x_cgs.x << " " << std::setw(15) << x_cgs.y << " " << std::setw(15)
           << x_cgs.z << "\n";

  for (int i = 0; i < STEPS_LIMIT; i++) {
    boris_step(t, x, v, q, m, E, B, dt, 1, x, v);
    t += dt;

    log_advancement(i, STEPS_LIMIT);

    // TODO: Fix for N particles
    x_cgs = space_to_cgs_units(x[0]);
    solution << std::setw(6) << time_to_cgs_units(t) << " " << std::setw(15)
             << x_cgs.x << " " << std::setw(15) << x_cgs.y << " "
             << std::setw(15) << x_cgs.z << "\n";
  }

  solution.close();

  std::cout << "Simulation complete!" << std::endl;

  return 0;
}
