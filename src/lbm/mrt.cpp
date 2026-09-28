#include "lbm/mrt.h"

namespace cglbm {
namespace lbm {

void mrt_rates(double s_nu, double s_e, double s_eps, double s_q, double* rates) {
    rates[kMomentRho] = 0.0;
    rates[kMomentEnergy] = s_e;
    rates[kMomentEnergySquare] = s_eps;
    rates[kMomentJx] = s_nu;
    rates[kMomentQx] = s_q;
    rates[kMomentJy] = s_nu;
    rates[kMomentQy] = s_q;
    rates[kMomentPxx] = s_nu;
    rates[kMomentPxy] = s_nu;
}

void mrt_collide(const double* f,
                 const double* f_eq,
                 const double* source,
                 const double* rates,
                 double dt,
                 double* out) {
    double change[kQ];
    for (int a = 0; a < kQ; a++) {
        double m_neq = 0.0, m_source = 0.0;
        for (int i = 0; i < kQ; i++) {
            m_neq += kMoment[a][i] * (f[i] - f_eq[i]);
            m_source += kMoment[a][i] * source[i];
        }
        change[a] = (-rates[a] * m_neq + (1.0 - 0.5 * rates[a]) * m_source * dt) / kMomentNorm[a];
    }
    for (int i = 0; i < kQ; i++) {
        double post = f[i];
        for (int a = 0; a < kQ; a++) {
            post += kMoment[a][i] * change[a];
        }
        out[i] = post;
    }
}

void add_third_moment_source(
    double dqx_dx, double dqy_dy, double s_e, double s_nu, double dt, double* out) {
    const double trace = 3.0 * (1.0 - 0.5 * s_e) * (dqx_dx + dqy_dy) * dt;
    const double normal = (1.0 - 0.5 * s_nu) * (dqx_dx - dqy_dy) * dt;
    for (int i = 0; i < kQ; i++) {
        out[i] += kMoment[kMomentEnergy][i] * trace / kMomentNorm[kMomentEnergy] +
                  kMoment[kMomentPxx][i] * normal / kMomentNorm[kMomentPxx];
    }
}

}  // namespace lbm
}  // namespace cglbm
