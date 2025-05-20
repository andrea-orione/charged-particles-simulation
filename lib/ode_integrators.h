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
 * @param y [in/out] An array containing the starting state.
 * It will also be used to store the final state.
 * @param dydt An array containing the derivative of the starting state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt The step of the increment.
 * @param nVars The number of equations in the system.
 */
void RK4_step(double t, double *y, double *dydt, Rhs_f rhs_f, double dt,
              int nVars);

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
 * @param y [in/out] An array containing the starting state.
 * It will also be used to store the final state.
 * @param dydt An array containing the derivative of the starting state.
 * @param rhs_f The right hand side of the differential equation.
 * @param dt The step of the increment.
 * @param nVars The number of equations in the system.
 * Maximum number supported is 256.
 *
 * @throws std::invalid_argument if nVars exceeds the supported value
 */
void adaptive_RK4_step(double t, double *y, Rhs_f rhs_f, int nVars);

#endif
