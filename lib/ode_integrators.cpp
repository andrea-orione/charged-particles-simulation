#include "ode_integrators.h"

#define MAX_NVARS 256

void RK4_step(const double t, double *y, double *dydt, Rhs_f rhs_f,
              const double dt, const int nVars) {
  double *yTemp = new double[nVars]; // this one gets reused
  double *k1 = new double[nVars];
  double *k2 = new double[nVars];
  double *k3 = new double[nVars];

  rhs_f(t, y, dydt);
  int i;
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + 0.5 * dt * dydt[i];

  rhs_f(t + 0.5 * dt, yTemp, k1);
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + 0.5 * dt * k1[i];

  rhs_f(t + 0.5 * dt, yTemp, k2);
  for (i = 0; i < nVars; i++)
    yTemp[i] = y[i] + dt * k2[i];

  rhs_f(t + dt, yTemp, k3);
  for (i = 0; i < nVars; i++)
    y[i] += dt * (dydt[i] + k3[i] + 2. * (k1[i] + k2[i])) / 6.;

  delete[] yTemp;
  delete[] k1;
  delete[] k2;
  delete[] k3;
}
