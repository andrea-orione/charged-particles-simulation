#ifndef ODE_INTEGRATORS
#define ODE_INTEGRATORS

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
 * midpoint of order 2 method.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) at a certain time, this returns the solution
 * y(t+dt) after a time-step dt.
 * This is done evaluating rhs twice.
 * The error goes with O(dt^3).
 *
 * @param t The time at which the step starts.
 * @param y [in/out] An array containing the starting state.
 * It will also be used to store the final state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt The step of the increment.
 * @param nVars The number of equations in the system.
 * Maximum number supported is 256.
 *
 * @throws std::invalid_argument if nVars exceeds the supported value
 */
void RK2_step(double t, double *y, Rhs_f rhs_f, double dt, int nVars);

/**
 * Calculate the solution of the ode after a time-step using Runge-Kutta
 * simpson of order 4 method.
 *
 * Given an ordinary differential equation of the form
 * dy/dt = rhs(y,t)
 * and given the solution y(t) at a certain time, this returns the solution
 * y(t+dt) after a time-step dt.
 * This is done evaluating rhs four times.
 * The error goes with O(dt^5).
 *
 * @param t The time at which the step starts.
 * @param y [in/out] An array containing the starting state.
 * It will also be used to store the final state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt The step of the increment.
 * @param nVars The number of equations in the system.
 * Maximum number supported is 256.
 *
 * @throws std::invalid_argument if nVars exceeds the supported value
 */
void RK4_step(double t, double *y, Rhs_f rhs_f, double dt, int nVars);

#endif
