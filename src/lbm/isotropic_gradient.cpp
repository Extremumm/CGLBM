#include "lbm/isotropic_gradient.h"

#include <cstring>

namespace cglbm {
namespace lbm {

namespace {

// The shell weights below already contain the 1/c_s^2 factor: they are
// normalised so that sum_l W c_a c_b = delta_ab, which makes the stencil exact
// for a linear field without any further scaling. For E4 that gives
// W(1) = 1/3 and W(2) = 1/12, which is w_k / c_s^2 with the D2Q9 weights
// w_k = 1/9 and 1/36 -- so E4 is exactly the stencil the solvers already used.

/// E4: shells |c|^2 = 1, 2. Isotropic to 4th order. [Sbragaglia 2007]
const StencilPoint kE4[] = {
    // |c|^2 = 1, W = 1/3
    {1, 0, 1.0 / 3.0},
    {-1, 0, 1.0 / 3.0},
    {0, 1, 1.0 / 3.0},
    {0, -1, 1.0 / 3.0},
    // |c|^2 = 2, W = 1/12
    {1, 1, 1.0 / 12.0},
    {-1, 1, 1.0 / 12.0},
    {1, -1, 1.0 / 12.0},
    {-1, -1, 1.0 / 12.0},
};

/// E6: shells |c|^2 = 1, 2, 4. Isotropic to 6th order. [Sbragaglia 2007, App. C]
const StencilPoint kE6[] = {
    // |c|^2 = 1, W = 4/15
    {1, 0, 4.0 / 15.0},
    {-1, 0, 4.0 / 15.0},
    {0, 1, 4.0 / 15.0},
    {0, -1, 4.0 / 15.0},
    // |c|^2 = 2, W = 1/10
    {1, 1, 1.0 / 10.0},
    {-1, 1, 1.0 / 10.0},
    {1, -1, 1.0 / 10.0},
    {-1, -1, 1.0 / 10.0},
    // |c|^2 = 4, W = 1/120
    {2, 0, 1.0 / 120.0},
    {-2, 0, 1.0 / 120.0},
    {0, 2, 1.0 / 120.0},
    {0, -2, 1.0 / 120.0},
};

/// E8: shells |c|^2 = 1, 2, 4, 5, 8. Isotropic to 8th order. [arXiv:2505.23647]
const StencilPoint kE8[] = {
    // |c|^2 = 1, W = 4/21
    {1, 0, 4.0 / 21.0},
    {-1, 0, 4.0 / 21.0},
    {0, 1, 4.0 / 21.0},
    {0, -1, 4.0 / 21.0},
    // |c|^2 = 2, W = 4/45
    {1, 1, 4.0 / 45.0},
    {-1, 1, 4.0 / 45.0},
    {1, -1, 4.0 / 45.0},
    {-1, -1, 4.0 / 45.0},
    // |c|^2 = 4, W = 1/60
    {2, 0, 1.0 / 60.0},
    {-2, 0, 1.0 / 60.0},
    {0, 2, 1.0 / 60.0},
    {0, -2, 1.0 / 60.0},
    // |c|^2 = 5, W = 2/315
    {1, 2, 2.0 / 315.0},
    {-1, 2, 2.0 / 315.0},
    {1, -2, 2.0 / 315.0},
    {-1, -2, 2.0 / 315.0},
    {2, 1, 2.0 / 315.0},
    {-2, 1, 2.0 / 315.0},
    {2, -1, 2.0 / 315.0},
    {-2, -1, 2.0 / 315.0},
    // |c|^2 = 8, W = 1/5040
    {2, 2, 1.0 / 5040.0},
    {-2, 2, 1.0 / 5040.0},
    {2, -2, 1.0 / 5040.0},
    {-2, -2, 1.0 / 5040.0},
};

int count_of(const StencilPoint* first, const StencilPoint* last) {
    return static_cast<int>(last - first);
}

}  // namespace

const StencilPoint* stencil_points(GradientStencil stencil, int* count) {
    switch (stencil) {
        case GradientStencil::E6:
            if (count != nullptr) {
                *count = count_of(kE6, kE6 + sizeof(kE6) / sizeof(kE6[0]));
            }
            return kE6;
        case GradientStencil::E8:
            if (count != nullptr) {
                *count = count_of(kE8, kE8 + sizeof(kE8) / sizeof(kE8[0]));
            }
            return kE8;
        case GradientStencil::E4:
        default:
            if (count != nullptr) {
                *count = count_of(kE4, kE4 + sizeof(kE4) / sizeof(kE4[0]));
            }
            return kE4;
    }
}

int stencil_reach(GradientStencil stencil) {
    return (stencil == GradientStencil::E4) ? 1 : 2;
}

const char* stencil_name(GradientStencil stencil) {
    switch (stencil) {
        case GradientStencil::E6:
            return "E6";
        case GradientStencil::E8:
            return "E8";
        case GradientStencil::E4:
        default:
            return "E4";
    }
}

bool stencil_from_name(const char* name, GradientStencil* stencil) {
    if (name == nullptr || stencil == nullptr) {
        return false;
    }
    if (std::strcmp(name, "E4") == 0 || std::strcmp(name, "e4") == 0) {
        *stencil = GradientStencil::E4;
        return true;
    }
    if (std::strcmp(name, "E6") == 0 || std::strcmp(name, "e6") == 0) {
        *stencil = GradientStencil::E6;
        return true;
    }
    if (std::strcmp(name, "E8") == 0 || std::strcmp(name, "e8") == 0) {
        *stencil = GradientStencil::E8;
        return true;
    }
    return false;
}

void gradient_periodic(const double* field,
                       int nx,
                       int ny,
                       int i,
                       int j,
                       GradientStencil stencil,
                       double* grad_x,
                       double* grad_y) {
    int count = 0;
    const StencilPoint* points = stencil_points(stencil, &count);

    double gx = 0.0;
    double gy = 0.0;
    for (int n = 0; n < count; ++n) {
        // (v % m + m) % m keeps the index non-negative for offsets of either
        // sign, and for a reach of 2 on a lattice at least 2 nodes wide.
        const int ip = ((i + points[n].cx) % nx + nx) % nx;
        const int jp = ((j + points[n].cy) % ny + ny) % ny;
        const double value = field[ip * ny + jp];
        gx += points[n].weight * points[n].cx * value;
        gy += points[n].weight * points[n].cy * value;
    }

    *grad_x = gx;
    *grad_y = gy;
}

void gradient_wall_y(const double* field,
                     int nx,
                     int ny,
                     int i,
                     int j,
                     GradientStencil stencil,
                     double* grad_x,
                     double* grad_y) {
    int count = 0;
    const StencilPoint* points = stencil_points(stencil, &count);

    double gx = 0.0;
    double gy = 0.0;
    for (int n = 0; n < count; ++n) {
        const int jp = j + points[n].cy;
        if (jp < 0 || jp >= ny) {
            // Past the wall there is no fluid node to read; the term is dropped
            // rather than wrapped, matching the solvers' own treatment.
            continue;
        }
        const int ip = ((i + points[n].cx) % nx + nx) % nx;
        const double value = field[ip * ny + jp];
        gx += points[n].weight * points[n].cx * value;
        gy += points[n].weight * points[n].cy * value;
    }

    *grad_x = gx;
    *grad_y = gy;
}

}  // namespace lbm
}  // namespace cglbm
