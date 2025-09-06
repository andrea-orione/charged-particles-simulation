#ifndef CONVERSIONS
#define CONVERSIONS

#include "init.h"
#include "vectors.h"

// TODO: Finish documentation
/**
 * This file contains the constants and functions to convert from physical units
 * to simulation units.
 *
 * The conversions are made based on the scale value set in `init.h`
 */

#ifndef SCALE
#define SCALE EARTH
#endif

constexpr double LIGHT_SPEED = 2.99792458e10;      // cm/s
constexpr double EARTH_MEAN_RADIUS = 6.3710088e8;  // cm
constexpr double ELECTRON_MASS = 9.1093837139e-26; // g
constexpr double ELECTRON_CHARGE = 4.8032047e-10;  // statC

// TODO: Add other scales
#if SCALE == EARTH
constexpr double unit_mag_field = 1.; // Gauss
#elif SCALE == LAB
constexpr double unit_mag_field = 1.e4; // Gauss = 1T
#endif

// TODO: Consider if is better to define another speed scale
// NOTE: Remember than in the integration algorithms it is not divided by c
constexpr double unit_mass = ELECTRON_MASS;     // g
constexpr double unit_charge = ELECTRON_CHARGE; // statC
constexpr double unit_speed = LIGHT_SPEED;      // cm/s

// TODO: Check if it works (not sure, doesn't feel right)
constexpr double unit_space = unit_speed * unit_speed * unit_mass /
                              (unit_charge * unit_mag_field); // Gauss
constexpr double unit_time = unit_space / unit_speed;         // s
constexpr double unit_el_field = unit_mag_field;              // statV/cm

constexpr double mass_to_sim_units(double mass) { return mass / unit_mass; }
constexpr double time_to_sim_units(double time) { return time / unit_time; }
constexpr double space_to_sim_units(double space) { return space / unit_space; }
constexpr double speed_to_sim_units(double speed) { return speed / unit_speed; }
constexpr double charge_to_sim_units(double charge) {
  return charge / unit_charge;
}

constexpr Vec3 space_to_sim_units(Vec3 space) {
  return scale(space, 1. / unit_space);
}
constexpr Vec3 speed_to_sim_units(Vec3 speed) {
  return scale(speed, 1. / unit_speed);
}
constexpr Vec3 mag_field_to_sim_units(Vec3 B) {
  return scale(B, 1. / unit_mag_field);
}
constexpr Vec3 el_field_to_sim_units(Vec3 E) {
  return scale(E, 1. / unit_el_field);
}

constexpr double mass_to_cgs_units(double mass) { return mass * unit_mass; }
constexpr double time_to_cgs_units(double time) { return time * unit_time; }
constexpr double space_to_cgs_units(double space) { return space * unit_space; }
constexpr double speed_to_cgs_units(double speed) { return speed * unit_speed; }
constexpr double charge_to_cgs_units(double charge) {
  return charge * unit_charge;
}

constexpr Vec3 space_to_cgs_units(Vec3 space) {
  return scale(space, unit_space);
}
constexpr Vec3 speed_to_cgs_units(Vec3 speed) {
  return scale(speed, unit_speed);
}
constexpr Vec3 mag_field_to_cgs_units(Vec3 B) {
  return scale(B, unit_mag_field);
}
constexpr Vec3 el_field_to_cgs_units(Vec3 E) { return scale(E, unit_el_field); }

#endif
