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
 * Calculate the solution of the ode after a time-step using Runge-Kutta
 * simpson of order 4 method.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) and its derivative dy/dt at a certain time,
 * this returns the solution y(t+dt) after a time-step dt.
 * This is done evaluating rhs three times.
 * The error goes with O(dt^5).
 * NOTE: The algorithm takes the derivative at the starting point so
 * if you need to run the step multiple times varying the timestep, the starting
 * derivative doesn't get recomputed needlessly.
 *
 * @param t The time at which the step starts.
 * @param y An array containing the starting state.
 * @param dydt An array containing the derivative of the starting state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt The step of the increment.
 * @param nVars The number of equations in the system.
 * @param y_out [out] An array where the next state will be saved.
 * It does not need to be distinct from y.
 */
void RK4_step(double t, double *y, double *dydt, Rhs_f rhs_f, double dt,
              int nVars, double *y_out);

/**
 * Calculate the solution of the ode after using Runge-Kutta simpson of order 4
 * method. The solution is calculated after a time-step such that the precision
 * of the solution is above a certain threshold or if the maximum number of
 * iteration is reached.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) at a certain time, this returns the solution
 * y(t+dt) after a time-step dt.
 * This process is repeated, each itme halfing the time-step, until the gain in
 * precision goes under a specified threshold or until a maximum number of
 * iteration is reached.
 * The gain in precision between two iterations is calculated as the sum of
 * of the absoluted value of the differences between the two iterations.
 * This is done evaluating rhs three times for each iteration.
 * The returned state is the one obtained from the last iteration, leading to a
 * possibli better solution.
 *
 * @param t The time at which the step starts.
 * @param y An array containing the starting state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt_init The initial time step of the
 * @param precis The minimum precision required.
 * @param max_iter The maximum number of iteration allowed
 * @param nVars The number of equations in the system.
 * @param y_out [out] An array where the next state will be saved.
 * It does not need to be distinct from y.
 *
 * @return The number of iteration done
 */
int adaptive_RK4_step(double t, double *y, Rhs_f rhs_f, double dt_init,
                      double precis, int max_iter, int nVars, double *y_out);

/**
 * Calculate the solution of the ode after using Runge-Kutta simpson of order 4
 * method. The solution is calculated after a time-step such that the precision
 * of the solution is above a certain threshold or if the maximum number of
 * iteration is reached.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) at a certain time, this returns the solution
 * y(t+dt) after a time-step dt.
 * This process is repeated, each itme halfing the time-step, until the gain in
 * precision goes under a specified threshold or until a maximum number of
 * iteration is reached.
 * The gain in precision between two iterations is calculated as the sum of
 * of the absoluted value of the differences between the two iterations.
 * This is done evaluating rhs three times for each iteration.
 * The returned state is the one obtained from the last iteration, leading to a
 * possibli better solution.
 *
 * @param t The time at which the step starts.
 * @param y An array containing the starting state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt_init The initial time step of the
 * @param precis The minimum precision required.
 * @param max_iter The maximum number of iteration allowed
 * @param nVars The number of equations in the system.
 * @param y_out [out] An array where the next state will be saved.
 * It does not need to be distinct from y.
 *
 * @return The number of iteration done
 */
void boris_step(Vec3 *x, Vec3 *v, double *q, Vec3Field E, Vec3Field B,
                double dt, int nParticles, Vec3 *x_out, Vec3 *y_out);

#endif
