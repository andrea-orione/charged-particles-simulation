#include "init.h"

#include "conversions.h"
#include "vectors.h"

// WARN: Always write in simulation units
void initialize_starting_state(Vec3 *start_pos, Vec3 *start_vel,
                               double *charges, double *masses) {
  for (int i=0; i< N_PART; i++) {
  start_pos[i] = space_to_sim_units(Vec3{30.-i, 0., 0.});
  start_vel[i] = speed_to_sim_units(Vec3{0., LIGHT_SPEED, 0.});
  charges[i] = charge_to_sim_units(ELECTRON_CHARGE);
  masses[i] = mass_to_sim_units(ELECTRON_MASS);
  }
}

Vec3 E(Vec3 x, double t) {
  double d = sqrt(x.x * x.x + x.y * x.y);
  d = 6.e2 / (d * d * d);
  return el_field_to_sim_units(Vec3{-x.x * d, -x.y * d, 0.});
}

Vec3 B(Vec3 x, double t) {
  return mag_field_to_sim_units(
      Vec3{0., 0., 1.e4 * sqrt(x.x * x.x + x.y * x.y)});
}
