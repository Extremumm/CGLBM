#include "lbm/quasi_static_mhd_3d.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cglbm {
namespace lbm {

namespace {

/// Nodes summed by one partial of a reduction.
///
/// The partition depends on the lattice size and on nothing else, so a dot
/// product sums the same numbers in the same order whatever the thread count.
/// Without that, every reduction below would make the run depend on how many
/// threads it happened to get, which is the one property the rest of this
/// solver goes out of its way to keep.
constexpr std::size_t kReductionChunk = 4096;

/// The D3Q7 lattice the potential march runs on: rest plus the six axes.
constexpr int kQ7 = 7;
constexpr int kXi7[kQ7][3] = {{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                              {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
constexpr int kOpposite7[kQ7] = {0, 2, 1, 4, 3, 6, 5};

/// Rest weight of D3Q7, which fixes `c_s^2 = (1 - w_0) / 3 = 1/4`.
constexpr double kW7Rest = 1. / 4.;
constexpr double kW7Axial = (1.0 - kW7Rest) / 6.0;
constexpr double kCs2Seven = (1.0 - kW7Rest) / 3.0;

/// The magic parameter of the two-relaxation-time collision.
///
/// `Lambda = (tau^+ - 1/2)(tau^- - 1/2)`, held at 1/4 so that where `sigma`
/// jumps the effective position of the jump does not move with the relaxation
/// time. With a single relaxation time it does, and here the relaxation time
/// spans the whole conductivity ratio across four nodes of interface.
constexpr double kMagicParameter = 1. / 4.;

constexpr double weight_seven(int q) {
    return q == 0 ? kW7Rest : kW7Axial;
}

}  // namespace

bool potential_solver_from_name(const char* name, PotentialSolver* solver) {
    if (std::strcmp(name, "fv") == 0 || std::strcmp(name, "finite-volume") == 0) {
        *solver = PotentialSolver::FiniteVolume;
        return true;
    }
    if (std::strcmp(name, "lbm") == 0 || std::strcmp(name, "lattice-boltzmann") == 0) {
        *solver = PotentialSolver::LatticeBoltzmann;
        return true;
    }
    return false;
}

const char* potential_solver_name(PotentialSolver solver) {
    return solver == PotentialSolver::LatticeBoltzmann ? "lbm" : "fv";
}

double magnetic_damping_time(double density, double conductivity, const double* b) {
    const double b_squared = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    const double rate = conductivity * b_squared;
    if (!(rate > 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    return density / rate;
}

QuasiStaticMhd3D::QuasiStaticMhd3D(
    int nx, int ny, int nz, bool wall_y, MhdPhysics physics, bool parallel)
    : nx_(nx), ny_(ny), nz_(nz), wall_y_(wall_y), parallel_(parallel),
      physics_(std::move(physics)) {
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0) {
        throw std::invalid_argument("lattice must have at least one node on each axis");
    }
    if (physics_.conductivity1 < 0.0 || physics_.conductivity2 < 0.0) {
        throw std::invalid_argument("both conductivities must be non-negative");
    }
    if (physics_.max_iterations <= 0) {
        throw std::invalid_argument("the potential solve needs at least one iteration");
    }

    // Nothing here pins the potential: every boundary is periodic or
    // insulating, so the operator keeps the constant nullspace and the solve
    // projects it out. A perfectly conducting wall would set this false.
    singular_ = true;

    nodes_ = static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) *
             static_cast<std::size_t>(nz_);
    chunk_count_ = (nodes_ + kReductionChunk - 1) / kReductionChunk;

    phi_ = Field3D(nx_, ny_, nz_);
    sigma_face_ = Field3D(nx_, ny_, nz_, 3);
    sigma_node_ = Field3D(nx_, ny_, nz_);
    if (physics_.solver == PotentialSolver::LatticeBoltzmann) {
        g_ = Field3D(nx_, ny_, nz_, kQ7);
        g_next_ = Field3D(nx_, ny_, nz_, kQ7);
        phi_previous_ = Field3D(nx_, ny_, nz_);
    }
    drive_ = Field3D(nx_, ny_, nz_, 3);
    emf_ = Field3D(nx_, ny_, nz_, 3);
    face_current_ = Field3D(nx_, ny_, nz_, 3);
    current_ = Field3D(nx_, ny_, nz_, 3);

    rhs_.assign(nodes_, 0.0);
    solution_.assign(nodes_, 0.0);
    residual_vector_.assign(nodes_, 0.0);
    direction_.assign(nodes_, 0.0);
    operator_direction_.assign(nodes_, 0.0);
    preconditioned_.assign(nodes_, 0.0);
    diagonal_.assign(nodes_, 1.0);
    partials_.assign(chunk_count_, 0.0);
}

namespace {

/// Node index into a flat `nx * ny * nz` array.
inline std::size_t node_index(int i, int j, int k, int ny, int nz) {
    return (static_cast<std::size_t>(i) * static_cast<std::size_t>(ny) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nz) +
           static_cast<std::size_t>(k);
}

/// The neighbour one step along `axis`, wrapped. A wall shows up as a face of
/// zero conductivity, not as a missing neighbour, so the wrap is unconditional.
inline void
step_along(int axis, int delta, int nx, int ny, int nz, int i, int j, int k, int* out) {
    out[0] = i;
    out[1] = j;
    out[2] = k;
    out[axis] += delta;
    const int extent[3] = {nx, ny, nz};
    out[axis] = (out[axis] + extent[axis]) % extent[axis];
}

}  // namespace

double QuasiStaticMhd3D::dot(const std::vector<double>& a, const std::vector<double>& b) const {
    const bool parallel = parallel_;
    const long chunks = static_cast<long>(chunk_count_);
    const std::size_t nodes = nodes_;
#pragma omp parallel for if (parallel)
    for (long c = 0; c < chunks; ++c) {
        const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
        const std::size_t end = std::min(begin + kReductionChunk, nodes);
        double sum = 0.0;
        for (std::size_t n = begin; n < end; ++n) {
            sum += a[n] * b[n];
        }
        partials_[static_cast<std::size_t>(c)] = sum;
    }
    double total = 0.0;
    for (std::size_t c = 0; c < chunk_count_; ++c) {
        total += partials_[c];
    }
    return total;
}

void QuasiStaticMhd3D::remove_mean(std::vector<double>& v) const {
    const bool parallel = parallel_;
    const long chunks = static_cast<long>(chunk_count_);
    const std::size_t nodes = nodes_;
#pragma omp parallel for if (parallel)
    for (long c = 0; c < chunks; ++c) {
        const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
        const std::size_t end = std::min(begin + kReductionChunk, nodes);
        double sum = 0.0;
        for (std::size_t n = begin; n < end; ++n) {
            sum += v[n];
        }
        partials_[static_cast<std::size_t>(c)] = sum;
    }
    double total = 0.0;
    for (std::size_t c = 0; c < chunk_count_; ++c) {
        total += partials_[c];
    }
    const double mean = total / static_cast<double>(nodes_);
#pragma omp parallel for if (parallel)
    for (long c = 0; c < chunks; ++c) {
        const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
        const std::size_t end = std::min(begin + kReductionChunk, nodes);
        for (std::size_t n = begin; n < end; ++n) {
            v[n] -= mean;
        }
    }
}

void QuasiStaticMhd3D::update_conductivity(const Field3D& phase) {
    const bool parallel = parallel_;
    const double sigma1 = physics_.conductivity1;
    const double sigma2 = physics_.conductivity2;
    const bool harmonic = physics_.harmonic_conductivity;
    const int nx = nx_, ny = ny_, nz = nz_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                // Volume fraction of component 1, which is what the phase field
                // already is once mapped from [-1, 1] onto [0, 1].
                const double here = std::min(1.0, std::max(0.0, 0.5 * (1.0 + phase(i, j, k))));
                // The node value the lattice-Boltzmann march relaxes at. The
                // finite-volume operator never reads it: there the conductivity
                // belongs on the face.
                sigma_node_(i, j, k) = here * sigma1 + (1.0 - here) * sigma2;
                for (int axis = 0; axis < 3; axis++) {
                    int n[3];
                    step_along(axis, 1, nx, ny, nz, i, j, k, n);
                    const double there =
                        std::min(1.0, std::max(0.0, 0.5 * (1.0 + phase(n[0], n[1], n[2]))));
                    const double fraction = 0.5 * (here + there);
                    double value;
                    if (harmonic) {
                        // Series conductance across the face, which is what a
                        // current crossing the interface sees.
                        const double denominator = fraction * sigma2 + (1.0 - fraction) * sigma1;
                        value = denominator > 0.0 ? sigma1 * sigma2 / denominator : 0.0;
                    } else {
                        value = fraction * sigma1 + (1.0 - fraction) * sigma2;
                    }
                    sigma_face_(i, j, k, axis) = value;
                }
            }
        }
    }
    if (wall_y_) {
        // The wrap-around face along y is the insulating wall. Dropping it from
        // the operator is the whole of `J . n = 0`, and it leaves the flux
        // balance of the two boundary layers of nodes intact.
#pragma omp parallel for collapse(2) if (parallel)
        for (int i = 0; i < nx_; i++) {
            for (int k = 0; k < nz_; k++) {
                sigma_face_(i, ny_ - 1, k, 1) = 0.0;
            }
        }
    }

    // The Jacobi preconditioner, and the row sum the operator below repeats.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double sum = 0.0;
                for (int axis = 0; axis < 3; axis++) {
                    int back[3];
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    sum += sigma_face_(i, j, k, axis);
                    sum += sigma_face_(back[0], back[1], back[2], axis);
                }
                // A node no current can reach leaves its row empty; the
                // preconditioner must not divide by that.
                diagonal_[node_index(i, j, k, ny, nz)] = sum > 0.0 ? sum : 1.0;
            }
        }
    }
}

void QuasiStaticMhd3D::update_drive(const Field3D& velocity) {
    const bool parallel = parallel_;
    const double bx = physics_.b[0], by = physics_.b[1], bz = physics_.b[2];
    const int nx = nx_, ny = ny_, nz = nz_;

    // The motional electromotive force u x B_0 at the nodes.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double ux = velocity(i, j, k, 0);
                const double uy = velocity(i, j, k, 1);
                const double uz = velocity(i, j, k, 2);
                emf_(i, j, k, 0) = uy * bz - uz * by;
                emf_(i, j, k, 1) = uz * bx - ux * bz;
                emf_(i, j, k, 2) = ux * by - uy * bx;
            }
        }
    }

    // sigma_f times the face-normal component, averaged onto the face.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                for (int axis = 0; axis < 3; axis++) {
                    int n[3];
                    step_along(axis, 1, nx, ny, nz, i, j, k, n);
                    const double face =
                        0.5 * (emf_(i, j, k, axis) + emf_(n[0], n[1], n[2], axis));
                    drive_(i, j, k, axis) = sigma_face_(i, j, k, axis) * face;
                }
            }
        }
    }

    // A phi = -div (sigma (u x B_0)); see the header for the sign.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double divergence = 0.0;
                for (int axis = 0; axis < 3; axis++) {
                    int back[3];
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    divergence += drive_(i, j, k, axis) - drive_(back[0], back[1], back[2], axis);
                }
                rhs_[node_index(i, j, k, ny, nz)] = -divergence;
            }
        }
    }
}

void QuasiStaticMhd3D::apply_operator(const Field3D& in, Field3D& out) const {
    const bool parallel = parallel_;
    const int nx = nx_, ny = ny_, nz = nz_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double here = in(i, j, k);
                double sum = 0.0;
                for (int axis = 0; axis < 3; axis++) {
                    int ahead[3], back[3];
                    step_along(axis, 1, nx, ny, nz, i, j, k, ahead);
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    sum += sigma_face_(i, j, k, axis) * (here - in(ahead[0], ahead[1], ahead[2]));
                    sum += sigma_face_(back[0], back[1], back[2], axis) *
                           (here - in(back[0], back[1], back[2]));
                }
                out(i, j, k) = sum;
            }
        }
    }
}

namespace {

/// `A v` on flat storage, sharing the stencil of `apply_operator`.
void apply_flat(const std::vector<double>& in,
                std::vector<double>& out,
                const Field3D& sigma_face,
                int nx,
                int ny,
                int nz,
                bool parallel) {
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                const double here = in[node_index(i, j, k, ny, nz)];
                double sum = 0.0;
                for (int axis = 0; axis < 3; axis++) {
                    int ahead[3], back[3];
                    step_along(axis, 1, nx, ny, nz, i, j, k, ahead);
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    sum += sigma_face(i, j, k, axis) *
                           (here - in[node_index(ahead[0], ahead[1], ahead[2], ny, nz)]);
                    sum += sigma_face(back[0], back[1], back[2], axis) *
                           (here - in[node_index(back[0], back[1], back[2], ny, nz)]);
                }
                out[node_index(i, j, k, ny, nz)] = sum;
            }
        }
    }
}

}  // namespace

void QuasiStaticMhd3D::solve_potential(const Field3D& rhs) {
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                rhs_[node_index(i, j, k, ny_, nz_)] = rhs(i, j, k);
            }
        }
    }
    solve_flat();
}

void QuasiStaticMhd3D::solve_flat() {
    const bool parallel = parallel_;
    const long chunks = static_cast<long>(chunk_count_);
    const std::size_t nodes = nodes_;
    std::vector<double>& z = preconditioned_;

    // The right-hand side is a discrete divergence, so it is already orthogonal
    // to the constant the operator cannot see. Removing what round-off left of
    // its mean keeps the whole iteration in that complement.
    if (singular_) {
        remove_mean(rhs_);
    }
    const double rhs_norm = std::sqrt(dot(rhs_, rhs_));
    if (!(rhs_norm > 0.0)) {
        // No drive: the potential is a constant, and the zero-mean
        // representative of it is zero.
        std::fill(solution_.begin(), solution_.end(), 0.0);
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    phi_(i, j, k) = 0.0;
                }
            }
        }
        iterations_ = 0;
        residual_ = 0.0;
        converged_ = true;
        return;
    }

    // The previous step's potential is the starting guess. The velocity moves
    // by a fraction of itself in one step, so after the first solve this is
    // worth an order of magnitude in iterations.
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                solution_[node_index(i, j, k, ny_, nz_)] = phi_(i, j, k);
            }
        }
    }
    if (singular_) {
        remove_mean(solution_);
    }

    apply_flat(solution_, operator_direction_, sigma_face_, nx_, ny_, nz_, parallel_);
#pragma omp parallel for if (parallel)
    for (long c = 0; c < chunks; ++c) {
        const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
        const std::size_t end = std::min(begin + kReductionChunk, nodes);
        for (std::size_t n = begin; n < end; ++n) {
            residual_vector_[n] = rhs_[n] - operator_direction_[n];
            z[n] = residual_vector_[n] / diagonal_[n];
        }
    }
    if (singular_) {
        remove_mean(residual_vector_);
        remove_mean(z);
    }
    direction_ = z;
    double rz = dot(residual_vector_, z);
    double residual_norm = std::sqrt(dot(residual_vector_, residual_vector_));

    const double threshold = physics_.tolerance * rhs_norm;
    int iteration = 0;
    while (residual_norm > threshold && iteration < physics_.max_iterations) {
        apply_flat(direction_, operator_direction_, sigma_face_, nx_, ny_, nz_, parallel_);
        if (singular_) {
            remove_mean(operator_direction_);
        }
        const double curvature = dot(direction_, operator_direction_);
        if (!(curvature > 0.0)) {
            break;  // the direction lies in the nullspace; nothing left to take
        }
        const double step = rz / curvature;
#pragma omp parallel for if (parallel)
        for (long c = 0; c < chunks; ++c) {
            const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
            const std::size_t end = std::min(begin + kReductionChunk, nodes);
            for (std::size_t n = begin; n < end; ++n) {
                solution_[n] += step * direction_[n];
                residual_vector_[n] -= step * operator_direction_[n];
                z[n] = residual_vector_[n] / diagonal_[n];
            }
        }
        if (singular_) {
            remove_mean(z);
        }
        const double rz_next = dot(residual_vector_, z);
        const double decay = rz_next / rz;
        rz = rz_next;
#pragma omp parallel for if (parallel)
        for (long c = 0; c < chunks; ++c) {
            const std::size_t begin = static_cast<std::size_t>(c) * kReductionChunk;
            const std::size_t end = std::min(begin + kReductionChunk, nodes);
            for (std::size_t n = begin; n < end; ++n) {
                direction_[n] = z[n] + decay * direction_[n];
            }
        }
        residual_norm = std::sqrt(dot(residual_vector_, residual_vector_));
        ++iteration;
    }

    if (singular_) {
        remove_mean(solution_);
    }
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                phi_(i, j, k) = solution_[node_index(i, j, k, ny_, nz_)];
            }
        }
    }
    iterations_ = iteration;
    residual_ = residual_norm / rhs_norm;
    converged_ = residual_norm <= threshold;
}

void QuasiStaticMhd3D::solve_lattice_boltzmann() {
    const bool parallel = parallel_;
    const bool wall = wall_y_;
    const int nx = nx_, ny = ny_, nz = nz_;

    // The same short-circuit the finite-volume path takes: with no drive the
    // potential is a constant, and its zero-mean representative is zero.
    double source_scale = 0.0;
    for (std::size_t n = 0; n < nodes_; ++n) {
        source_scale = std::max(source_scale, std::fabs(rhs_[n]));
    }
    if (!(source_scale > 0.0)) {
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    phi_(i, j, k) = 0.0;
                    for (int q = 0; q < kQ7; q++) {
                        g_(i, j, k, q) = 0.0;
                    }
                }
            }
        }
        iterations_ = 0;
        residual_ = 0.0;
        converged_ = true;
        return;
    }

    // Only the steady state is wanted, and it survives scaling the
    // conductivity and the source by the same constant. Scale so the largest
    // relaxation time is the one the case asked for, which is as fast as the
    // march can be run without leaving that fluid under-resolved.
    double sigma_max = 0.0;
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                sigma_max = std::max(sigma_max, sigma_node_(i, j, k));
            }
        }
    }
    if (!(sigma_max > 0.0)) {
        // Nothing conducts: no current, and the potential is immaterial.
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    phi_(i, j, k) = 0.0;
                }
            }
        }
        iterations_ = 0;
        residual_ = 0.0;
        converged_ = true;
        return;
    }
    const double tau_max = std::max(physics_.lbm_tau_max, 0.5 + 1e-6);
    const double scale = kCs2Seven * (tau_max - 0.5) / sigma_max;

    // Warm start: the populations are left from the previous step, and the
    // potential they carry is the previous answer. On the first call they are
    // zero, which is the zero potential.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double sum = 0.0;
                for (int q = 0; q < kQ7; q++) {
                    sum += g_(i, j, k, q);
                }
                phi_(i, j, k) = sum;
                phi_previous_(i, j, k) = sum;
            }
        }
    }

    const int check_every = std::max(1, physics_.lbm_check_interval);
    int sweep = 0;
    double change = 0.0;
    bool converged = false;
    while (sweep < physics_.max_iterations) {
        // Collide, two relaxation times, and stream in the same pass.
#pragma omp parallel for collapse(3) if (parallel)
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    const double diffusivity = scale * sigma_node_(i, j, k);
                    const double tau_minus = diffusivity / kCs2Seven + 0.5;
                    const double omega_minus = 1.0 / tau_minus;
                    // Lambda fixed, so tau^+ follows. A perfect insulator has
                    // tau^- = 1/2 and omega^+ = 0, which leaves the collision
                    // as a plain reversal: no flux, which is what it should be.
                    const double slack = tau_minus - 0.5;
                    const double omega_plus = slack > 0.0 ? 1.0 / (kMagicParameter / slack + 0.5)
                                                          : 0.0;
                    const double potential = phi_(i, j, k);
                    const double source = scale * rhs_[node_index(i, j, k, ny, nz)];

                    for (int q = 0; q < kQ7; q++) {
                        const int back = kOpposite7[q];
                        const double here = g_(i, j, k, q);
                        const double mirrored = g_(i, j, k, back);
                        const double symmetric = 0.5 * (here + mirrored);
                        const double antisymmetric = 0.5 * (here - mirrored);
                        const double weight = weight_seven(q);
                        const double post = here - omega_plus * (symmetric - weight * potential) -
                                            omega_minus * antisymmetric + weight * source;

                        int ip = (i + kXi7[q][0] + nx) % nx;
                        int kp = (k + kXi7[q][2] + nz) % nz;
                        int jp = j + kXi7[q][1];
                        int qp = q;
                        if (wall && (jp < 0 || jp >= ny)) {
                            // Half-way bounce-back: an insulating wall is a
                            // wall no population crosses, which is `J . n = 0`
                            // without a ghost node, exactly as the flow solver
                            // closes its own no-slip wall.
                            jp = j;
                            qp = kOpposite7[q];
                        } else if (!wall) {
                            jp = (j + kXi7[q][1] + ny) % ny;
                        }
                        g_next_(ip, jp, kp, qp) = post;
                    }
                }
            }
        }
        std::swap(g_, g_next_);
        ++sweep;

#pragma omp parallel for collapse(3) if (parallel)
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    double sum = 0.0;
                    for (int q = 0; q < kQ7; q++) {
                        sum += g_(i, j, k, q);
                    }
                    phi_(i, j, k) = sum;
                }
            }
        }

        if (sweep % check_every == 0) {
            double largest_change = 0.0;
            double largest = 0.0;
            for (int i = 0; i < nx_; i++) {
                for (int j = 0; j < ny_; j++) {
                    for (int k = 0; k < nz_; k++) {
                        const double value = phi_(i, j, k);
                        largest_change =
                            std::max(largest_change, std::fabs(value - phi_previous_(i, j, k)));
                        largest = std::max(largest, std::fabs(value));
                        phi_previous_(i, j, k) = value;
                    }
                }
            }
            change = largest > 0.0 ? largest_change / largest : largest_change;
            if (change <= physics_.tolerance) {
                converged = true;
                break;
            }
        }
    }

    if (singular_) {
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    solution_[node_index(i, j, k, ny_, nz_)] = phi_(i, j, k);
                }
            }
        }
        remove_mean(solution_);
        for (int i = 0; i < nx_; i++) {
            for (int j = 0; j < ny_; j++) {
                for (int k = 0; k < nz_; k++) {
                    phi_(i, j, k) = solution_[node_index(i, j, k, ny_, nz_)];
                }
            }
        }
    }
    iterations_ = sweep;
    residual_ = change;
    converged_ = converged;
}

void QuasiStaticMhd3D::update_current() {
    const bool parallel = parallel_;
    const int nx = nx_, ny = ny_, nz = nz_;

    // The face current, built on the difference the potential was solved
    // against. `sum_faces` of it is the residual of that solve and nothing
    // else, which is the property the whole arrangement exists for.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double here = phi_(i, j, k);
                for (int axis = 0; axis < 3; axis++) {
                    int ahead[3];
                    step_along(axis, 1, nx, ny, nz, i, j, k, ahead);
                    const double gradient = phi_(ahead[0], ahead[1], ahead[2]) - here;
                    face_current_(i, j, k, axis) =
                        drive_(i, j, k, axis) - sigma_face_(i, j, k, axis) * gradient;
                }
            }
        }
    }

    // Averaged back to the node, which is where the force is wanted.
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                for (int axis = 0; axis < 3; axis++) {
                    int back[3];
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    current_(i, j, k, axis) =
                        0.5 * (face_current_(i, j, k, axis) +
                               face_current_(back[0], back[1], back[2], axis));
                }
            }
        }
    }
}

void QuasiStaticMhd3D::solve(const Field3D& velocity, const Field3D& phase) {
    update_conductivity(phase);
    update_drive(velocity);
    if (physics_.solver == PotentialSolver::LatticeBoltzmann) {
        solve_lattice_boltzmann();
    } else {
        solve_flat();
    }
    update_current();
}

double QuasiStaticMhd3D::charge_imbalance() const {
    const int nx = nx_, ny = ny_, nz = nz_;
    double worst = 0.0;
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double divergence = 0.0;
                for (int axis = 0; axis < 3; axis++) {
                    int back[3];
                    step_along(axis, -1, nx, ny, nz, i, j, k, back);
                    divergence += face_current_(i, j, k, axis) -
                                  face_current_(back[0], back[1], back[2], axis);
                }
                worst = std::max(worst, std::fabs(divergence));
            }
        }
    }
    return worst;
}

void QuasiStaticMhd3D::add_lorentz_force(Field3D& force) const {
    const bool parallel = parallel_;
    const double bx = physics_.b[0], by = physics_.b[1], bz = physics_.b[2];
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double jx = current_(i, j, k, 0);
                const double jy = current_(i, j, k, 1);
                const double jz = current_(i, j, k, 2);
                force(i, j, k, 0) += jy * bz - jz * by;
                force(i, j, k, 1) += jz * bx - jx * bz;
                force(i, j, k, 2) += jx * by - jy * bx;
            }
        }
    }
}

}  // namespace lbm
}  // namespace cglbm
