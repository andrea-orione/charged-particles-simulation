#include "ode_integrators.h"
#include <cstring>

void RK4_step(const double t, double *y, double *dydt, Rhs_f rhs_f,
              const double dt, const int nVars, double *y_out) {
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
    y_out[i] = y[i] + dt * (dydt[i] + k3[i] + 2. * (k1[i] + k2[i])) / 6.;

  delete[] yTemp;
  delete[] k1;
  delete[] k2;
  delete[] k3;
}

int adaptive_RK4_step(const double t, double *y, Rhs_f rhs_f,
                      const double dt_init, const double precis,
                      const int max_iter, const int nVars, double *y_out) {
  int iter = 1;
  double *dydt_init = new double[nVars];
  double *y_next = new double[nVars];
  double *y_next_bet = new double[nVars];
  rhs_f(t, y, dydt_init);
  double dt = dt_init;

  RK4_step(t, y, dydt_init, rhs_f, dt, nVars, y_next_bet);

  double err = precis + 1.;
  while (err > precis && iter < max_iter) {
    memcpy(y_next, y_next_bet, nVars);

    iter++;
  }

  memcpy(y_out, y_next_bet, nVars);
  delete[] dydt_init;
  delete[] y_next;
  delete[] y_next_bet;

  return iter;
}

void boris_step(Vec3 *x, Vec3 *v, double q, Vec3Field E, Vec3Field B, double dt,
                int nParticles, Vec3 *x_out, Vec3 *y_out) {

  Vec3 *acc = new Vec3[nParticles];

  int i;
  for (i = 0; i < nVars; i++)
    x[i] += 0.5 * dt * v[i];

  accFunc(x, acc, nVars);
  for (i = 0; i < nVars; i++)
    v[i] = v[i] + dt * acc[i];

  for (i = 0; i < nVars; i++)
    x[i] += +0.5 * dt * v[i];
}
}
