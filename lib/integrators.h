#ifndef ODE_INTEGRATORS
#define ODE_INTEGRATORS

#include "vectors.h"

/**
 * The right hand side of a differential equation of the form
 * dy/dt = rhs(y,t)
 * It should have a prodotype rhs_f(double t, double *y, double *result) where
 *  - t is the time of evaluation
 *  - y is the solution at time t
 *  - result is an array where the right hand side will be stored
 */
typedef void (*Rhs_f)(double, double *, double *);

/**
 * Calculate the solution of a system of particles in an electro-magnetic field
 * after a time-step using Runge-Kutta simpson of order 4 method.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) and its derivative dy/dt at a certain time,
 * this returns the solution y(t+dt) after a time-step dt.
 * This is done evaluating rhs three times.
 * The error goes with O(dt^5).
 * NOTE: This function is specific for the type of system described above and
 * should not be used to solve other ODEs.
 *
 * @param t The time at which the step starts.
 * @param x The position of the particles at time t.
 * @param v The velocity of the particles.
 * @param q The charge of the particles.
 * @param m The mass of the particles.
 * @param E The electric field of the system.
 * @param B The acceleration of the system.
 * @param dt The step of the increment.
 * @param n_part The number of particles in the system
 * @param x_out [out] An array where the next positions will be saved.
 * It does not need to be distinct from y.
 * @param v_out [out] An array where the next positions will be saved.
 * It does not need to be distinct from y.
 */
void RK4_step(double t, Vec3 *x, Vec3 *v, double *q, double *m, Vec3Field E,
              Vec3Field B, double dt, int n_part, Vec3 *x_out, Vec3 *v_out);

/**
 * Calculates the state of a system of particles in an electro-magnetic field
 * after a time-step using Boris algorithm.
 *
 * Given a system with a state with 3D positions and 3D velocities of n
 * particles that evolve in an electro-magnetic field
 * dx/dt = v
 * dv/dt = q/m * (E + v/c x B)
 * and given the state of the system [x(t),v(t - dt/2)] this returns
 * the state after a time-step dt.
 * This is done evaluating E and B once each.
 * The error goes with O(dt^2).
 * This algorithm has the property of conserving the phase space volume and
 * therefore the energy of the system.
 * NOTE: For the conservation property to hold it is foundamental to evaluate
 * the position and the velocity in a staggered way.
 *
 * @param t The time at which the step starts.
 * @param x The position of the particles at time t.
 * @param v The velocity of the particles.
 * @param q The charge of the particles.
 * @param m The mass of the particles.
 * @param E The electric field of the system.
 * @param B The acceleration of the system.
 * @param dt The step of the increment.
 * @param n_part The number of particles in the system
 * @param x_out [out] An array where the next positions will be saved.
 * It does not need to be distinct from y.
 * @param v_out [out] An array where the next velocities will be saved.
 * It does not need to be distinct from y.
 */
void boris_step(double t, Vec3 *x, Vec3 *v, double *q, double *m, Vec3Field E,
                Vec3Field B, double dt, int n_part, Vec3 *x_out, Vec3 *v_out);

#endif
