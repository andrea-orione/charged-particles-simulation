#ifndef INIT
#define INIT

#include "vectors.h"

/**
 * This is the file for configuring the simulation.
 *
 * Every parameter of the simulation will be configured here.
 */

// WARN: Do not touch these lines
#define TIME_LIM_MODE 0 // If the simulation should stop after a time limit
#define STEP_LIM_MODE 1 // If the simulation should stop after n iterations

#define MODE STEP_LIM_MODE

#if MODE == TIME_LIM_MODE
/// The time after which the simulation should stop
/// TODO: Make considerations on units of measumemenst
constexpr double TIME_LIMIT = 90;
#elif MODE == STEP_LIM_MODE
/// The number of steps after the simulation should stop
constexpr int STEPS_LIMIT = 100000;
#endif

/// The time between a step and another in seconds
constexpr double STEP_SIZE = 0.03e-10;

// WARN: Do not touch these lines
#define EARTH 0
#define LAB 1
// TODO: Make units docs
// TODO: Add more scales
/**
 * This variable is used to set the scale used to convert into simulation units.
 * For more information, check out `units_docs.txt`
 * Possible values are:
 * - EARTH
 * Default is EARTH
 */
#define SCALE LAB

constexpr int N_PART = 4;
void initialize_starting_state(Vec3 *start_pos, Vec3 *start_vel,
                               double *charges, double *masses);
Vec3 E(Vec3 x, double t);
Vec3 B(Vec3 x, double t);

#endif
