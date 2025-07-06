#include "integrators.h"
#include "vectors.h"

void RK4_step(const double t, const Vec3 *const x, const Vec3 *const v,
              const double *const q, const double *const m, Vec3Field E,
              Vec3Field B, const double dt, const int n_part, Vec3 *const x_out,
              Vec3 *const v_out) {
  Vec3 x_tem, v_tem; // these ones gets reused
  Vec3 kx, kv;       // these one gets reused
  Vec3 kx_tot, kv_tot;

  double qm;
  double dt_2 = dt * 0.5;
  double dt_6 = dt / 6.;
  for (int i = 0; i < n_part; i++) {
    qm = q[i] / m[i];

    // K_1
    kx_tot = v[i];
    kv_tot = scale(sum(E(x[i], t), cross(v[i], B(x[i], t))), qm);

    // K_2
    x_tem = sum(x[i], scale(kx_tot, dt_2));
    v_tem = sum(v[i], scale(kv_tot, dt_2));

    kx = v_tem;
    kv = scale(sum(E(x_tem, t + dt_2), cross(v_tem, B(x_tem, t + dt_2))), qm);

    kx_tot = sum(kx_tot, scale(kx, 2));
    kv_tot = sum(kv_tot, scale(kv, 2));

    // K_3
    x_tem = sum(x[i], scale(kx, dt_2));
    v_tem = sum(v[i], scale(kv, dt_2));

    kx = v_tem;
    kv = scale(sum(E(x_tem, t + dt_2), cross(v_tem, B(x_tem, t + dt_2))), qm);

    kx_tot = sum(kx_tot, scale(kx, 2));
    kv_tot = sum(kv_tot, scale(kv, 2));

    // K_4
    x_tem = sum(x[i], scale(kx, dt));
    v_tem = sum(v[i], scale(kv, dt));

    kx = v_tem;
    kv = scale(sum(E(x_tem, t + dt), cross(v_tem, B(x_tem, t + dt))), qm);

    kx_tot = sum(kx_tot, kx);
    kv_tot = sum(kv_tot, kv);

    x_out[i] = sum(x[i], scale(kx_tot, dt_6));
    v_out[i] = sum(v[i], scale(kv_tot, dt_6));
  }
}

void boris_step(const double t, const Vec3 *const x, const Vec3 *const v,
                const double *const q, const double *const m, Vec3Field E,
                Vec3Field B, const double dt, const int n_part,
                Vec3 *const x_out, Vec3 *const v_out) {
  for (int i = 0; i < n_part; i++) {
    double qp = dt * q[i] / (2 * m[i]);
    Vec3 h = scale(B(x[i], t), qp);
    Vec3 qE = scale(E(x[i], t), qp);
    double h_mod2 = squared_mod(h);
    Vec3 s = scale(h, 2 / (1 + h_mod2));
    Vec3 u = sum(v[i], qE);
    Vec3 up = sum(u, cross(sum(u, cross(u, h)), s));

    v_out[i] = sum(up, qE);
    x_out[i] = sum(x[i], scale(v_out[i], dt));
  }
}
