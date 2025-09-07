#include "utils.h"

#include "conversions.h"
#include "vectors.h"

#include <cstdio>
#include <iomanip>
#include <iostream>

void log_units() {
  printf("---------------- Simulation units ----------------\n");
  printf("space_unit            = %.4e cm\n", unit_space);
  printf("speed_unit  [derived] = %.4e cm/s\n", unit_speed);
  printf("time_unit   [derived] = %.4e s\n\n", unit_time);

  printf("mass_unit             = %.4e g\n", unit_mass);
  printf("charge_unit           = %.4e statC\n\n", unit_charge);

  printf("B_unit                = %.4e G\n", unit_mag_field);
  printf("E_unit      [derived] = %.4e statV/cm\n", unit_el_field);
  printf("--------------------------------------------------\n");
}

void log_progress(int current_step, int tot_steps) {
  if (current_step == -1)
    current_step = 0;

  printf("Step %i of %i (%.2f%%)\r", current_step, tot_steps,
         100. * current_step / tot_steps);
  fflush(stdout);
}

void save_pos_state(std::ofstream *file, double t, Vec3 *xs, int n_part) {
  std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(7)
            << std::right;

  Vec3 x_cgs;
  *file << time_to_cgs_units(t) << ",";
  for (int i = 0; i < n_part-1; i++) {
    x_cgs = space_to_cgs_units(xs[i]);
    *file << x_cgs.x << "," << x_cgs.y << "," << x_cgs.z << ",";
  }
  x_cgs = space_to_cgs_units(xs[n_part-1]);
  *file << x_cgs.x << "," << x_cgs.y << "," << x_cgs.z << "\n";
}
