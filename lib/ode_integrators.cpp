#include "ode_integrators.h"

#include <stdexcept>

#define MAX_NVARS 256

void RK2_step(const double t, double *const y, Rhs_f rhs_f, const double dt,
              const int nVars) {
  // Preliminary check
  if (nVars > MAX_NVARS)
    throw std::invalid_argument(
        "rungeKutta2Step: nVars too high. Maximum supported is 256");

  // Reuse the k to save space
  double y1[MAX_NVARS], k[MAX_NVARS];

  rhs_f(t, y, k);
  int i;
  for (i = 0; i < nVars; i++)
    y1[i] = y[i] + 0.5 * dt * k[i];

  rhs_f(t + 0.5 * dt, y1, k);
  for (i = 0; i < nVars; i++)
    y[i] += dt * k[i];
}

void RK4_step(const double t, double *y, Rhs_f rhs_f, const double dt,
              const int nVars) {
  // Preliminary check
  if (nVars > MAX_NVARS)
    throw std::invalid_argument(
        "rungeKutta4Step: nVars too high. Maximum supported is 256");

  // Reuse the ytemp to save space
  double yTemp[MAX_NVARS];
  double k1[MAX_NVARS], k2[MAX_NVARS], k3[MAX_NVARS], k4[MAX_NVARS];

  rhs_f(t, y, k1);
  int i;
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + 0.5 * dt * k1[i];

  rhs_f(t + 0.5 * dt, yTemp, k2);
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + 0.5 * dt * k2[i];

  rhs_f(t + 0.5 * dt, yTemp, k3);
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + dt * k3[i];

  rhs_f(t + dt, yTemp, k4);
  for (i = 0; i < nVars; i++)
    y[i] += dt * (k1[i] + k4[i] + 2. * (k2[i] + k3[i])) / 6.;
}
