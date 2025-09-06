#include "init.h"

#include "conversions.h"
#include "vectors.h"

// WARN: Always write in simulation units
void initialize_starting_state(Vec3 *start_pos, Vec3 *start_vel,
                               double *charges, double *masses) {
  start_pos[0] = space_to_sim_units(Vec3{0., 0., 0.});
  start_vel[0] = speed_to_sim_units(Vec3{0., 0.01*LIGHT_SPEED, 0.});
  charges[0] = charge_to_sim_units(ELECTRON_CHARGE);
  masses[0] = mass_to_sim_units(ELECTRON_MASS);
}

Vec3 B(Vec3 x, double t) { return mag_field_to_sim_units(Vec3{0., 0., 1.}); }
Vec3 E(Vec3 x, double t) { return el_field_to_sim_units(Vec3{0., 0.5, 0.}); }
