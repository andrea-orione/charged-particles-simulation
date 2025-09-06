#include "utils.h"

#include "conversions.h"

#include <cstdio>

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

void log_advancement(int current_step, int tot_steps) {
  if (current_step == -1) current_step = 0;
  
  printf("Step %i of %i (%.2f%%)\r", current_step, tot_steps,
         100. * current_step / tot_steps);
  fflush(stdout);
}
