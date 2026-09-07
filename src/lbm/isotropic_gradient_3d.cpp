#include "lbm/isotropic_gradient_3d.h"

#include <cstring>

namespace cglbm {
namespace lbm {

namespace {

/// Build a shell of offsets with a given squared length, at a given weight.
///
/// Written out rather than generated so that the tables below are visible in
/// the source and can be read against the derivation in the header.
constexpr StencilPoint3D kE4[] = {
    // |c|^2 = 1, weight 1/6
    {1, 0, 0, 1. / 6.},
    {-1, 0, 0, 1. / 6.},
    {0, 1, 0, 1. / 6.},
    {0, -1, 0, 1. / 6.},
    {0, 0, 1, 1. / 6.},
    {0, 0, -1, 1. / 6.},
    // |c|^2 = 2, weight 1/12
    {1, 1, 0, 1. / 12.},
    {1, -1, 0, 1. / 12.},
    {-1, 1, 0, 1. / 12.},
    {-1, -1, 0, 1. / 12.},
    {1, 0, 1, 1. / 12.},
    {1, 0, -1, 1. / 12.},
    {-1, 0, 1, 1. / 12.},
    {-1, 0, -1, 1. / 12.},
    {0, 1, 1, 1. / 12.},
    {0, 1, -1, 1. / 12.},
    {0, -1, 1, 1. / 12.},
    {0, -1, -1, 1. / 12.},
};

constexpr StencilPoint3D kE6[] = {
    // |c|^2 = 1, weight 2/15
    {1, 0, 0, 2. / 15.},
    {-1, 0, 0, 2. / 15.},
    {0, 1, 0, 2. / 15.},
    {0, -1, 0, 2. / 15.},
    {0, 0, 1, 2. / 15.},
    {0, 0, -1, 2. / 15.},
    // |c|^2 = 2, weight 1/15
    {1, 1, 0, 1. / 15.},
    {1, -1, 0, 1. / 15.},
    {-1, 1, 0, 1. / 15.},
    {-1, -1, 0, 1. / 15.},
    {1, 0, 1, 1. / 15.},
    {1, 0, -1, 1. / 15.},
    {-1, 0, 1, 1. / 15.},
    {-1, 0, -1, 1. / 15.},
    {0, 1, 1, 1. / 15.},
    {0, 1, -1, 1. / 15.},
    {0, -1, 1, 1. / 15.},
    {0, -1, -1, 1. / 15.},
    // |c|^2 = 3, weight 1/60
    {1, 1, 1, 1. / 60.},
    {1, 1, -1, 1. / 60.},
    {1, -1, 1, 1. / 60.},
    {1, -1, -1, 1. / 60.},
    {-1, 1, 1, 1. / 60.},
    {-1, 1, -1, 1. / 60.},
    {-1, -1, 1, 1. / 60.},
    {-1, -1, -1, 1. / 60.},
    // |c|^2 = 4, weight 1/120
    {2, 0, 0, 1. / 120.},
    {-2, 0, 0, 1. / 120.},
    {0, 2, 0, 1. / 120.},
    {0, -2, 0, 1. / 120.},
    {0, 0, 2, 1. / 120.},
    {0, 0, -2, 1. / 120.},
};

}  // namespace

const StencilPoint3D* stencil_points_3d(GradientStencil3D stencil, int* count) {
    if (stencil == GradientStencil3D::E6) {
        *count = static_cast<int>(sizeof(kE6) / sizeof(kE6[0]));
        return kE6;
    }
    *count = static_cast<int>(sizeof(kE4) / sizeof(kE4[0]));
    return kE4;
}

int stencil_reach_3d(GradientStencil3D stencil) {
    return stencil == GradientStencil3D::E6 ? 2 : 1;
}

const char* stencil_name_3d(GradientStencil3D stencil) {
    return stencil == GradientStencil3D::E6 ? "E6" : "E4";
}

bool stencil_from_name_3d(const char* name, GradientStencil3D* stencil) {
    if (name == nullptr) {
        return false;
    }
    if (std::strcmp(name, "E4") == 0 || std::strcmp(name, "e4") == 0) {
        *stencil = GradientStencil3D::E4;
        return true;
    }
    if (std::strcmp(name, "E6") == 0 || std::strcmp(name, "e6") == 0) {
        *stencil = GradientStencil3D::E6;
        return true;
    }
    return false;
}

void gradient_periodic_3d(const double* field,
                          int nx,
                          int ny,
                          int nz,
                          int i,
                          int j,
                          int k,
                          GradientStencil3D stencil,
                          double* grad_x,
                          double* grad_y,
                          double* grad_z) {
    int count = 0;
    const StencilPoint3D* points = stencil_points_3d(stencil, &count);
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int p = 0; p < count; ++p) {
        const int ii = (i + points[p].cx % nx + nx) % nx;
        const int jj = (j + points[p].cy % ny + ny) % ny;
        const int kk = (k + points[p].cz % nz + nz) % nz;
        const double value = field[(static_cast<long>(ii) * ny + jj) * nz + kk];
        gx += points[p].weight * points[p].cx * value;
        gy += points[p].weight * points[p].cy * value;
        gz += points[p].weight * points[p].cz * value;
    }
    *grad_x = gx;
    *grad_y = gy;
    *grad_z = gz;
}

void gradient_wall_y_3d(const double* field,
                        int nx,
                        int ny,
                        int nz,
                        int i,
                        int j,
                        int k,
                        GradientStencil3D stencil,
                        double* grad_x,
                        double* grad_y,
                        double* grad_z) {
    int count = 0;
    const StencilPoint3D* points = stencil_points_3d(stencil, &count);
    double gx = 0.0, gy = 0.0, gz = 0.0;
    for (int p = 0; p < count; ++p) {
        const int jj = j + points[p].cy;
        if (jj < 0 || jj >= ny) {
            continue;  // beyond a wall: contributes nothing, as in two dimensions
        }
        const int ii = (i + points[p].cx % nx + nx) % nx;
        const int kk = (k + points[p].cz % nz + nz) % nz;
        const double value = field[(static_cast<long>(ii) * ny + jj) * nz + kk];
        gx += points[p].weight * points[p].cx * value;
        gy += points[p].weight * points[p].cy * value;
        gz += points[p].weight * points[p].cz * value;
    }
    *grad_x = gx;
    *grad_y = gy;
    *grad_z = gz;
}

}  // namespace lbm
}  // namespace cglbm
