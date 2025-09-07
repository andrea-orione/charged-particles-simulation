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

  std::ofstream solution;
  // TODO: Organize file output
  solution.open("../output/main.csv");

  log_units();

  double t = 0.;
  log_progress(-1, STEPS_LIMIT);
  save_pos_state(&solution, t, x, N_PART);

  for (int i = 0; i < STEPS_LIMIT; i++) {
    boris_step(t, x, v, q, m, E, B, dt, N_PART, x, v);
    t += dt;

    log_progress(i, STEPS_LIMIT);
    save_pos_state(&solution, t, x, N_PART);
  }

  solution.close();

  std::cout << "Simulation complete!" << std::endl;

  return 0;
}
