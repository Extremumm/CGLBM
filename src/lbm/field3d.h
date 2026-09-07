#ifndef CGLBM_LBM_FIELD3D_H
#define CGLBM_LBM_FIELD3D_H

#include <cstddef>
#include <vector>

/// Lattice-shaped storage in three dimensions.
///
/// The same object as `Field`, one axis wider. It is a separate class rather
/// than an extra dimension on `Field` because the two-dimensional solvers index
/// `Field` as `(i, j, k)` with `k` the velocity, and adding a `z` there would
/// have made every existing call site ambiguous.
///
/// Element `(i, j, k, q)` sits at `(((i * ny + j) * nz + k) * depth + q)`, so a
/// `depth == 1` field is a contiguous scalar lattice, which is what the
/// gradient stencils take. Storage is value-initialised to zero.

namespace cglbm {
namespace lbm {

class Field3D {
public:
    Field3D() = default;

    Field3D(int nx, int ny, int nz, int depth = 1)
        : nx_(nx), ny_(ny), nz_(nz), depth_(depth),
          data_(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                    static_cast<std::size_t>(nz) * static_cast<std::size_t>(depth),
                0.0) {}

    /// Value at node (i, j, k) of a scalar field.
    double& operator()(int i, int j, int k) {
        return data_[index(i, j, k, 0)];
    }
    double operator()(int i, int j, int k) const {
        return data_[index(i, j, k, 0)];
    }

    /// Component q at node (i, j, k).
    double& operator()(int i, int j, int k, int q) {
        return data_[index(i, j, k, q)];
    }
    double operator()(int i, int j, int k, int q) const {
        return data_[index(i, j, k, q)];
    }

    double* data() {
        return data_.data();
    }
    const double* data() const {
        return data_.data();
    }

    int nx() const {
        return nx_;
    }
    int ny() const {
        return ny_;
    }
    int nz() const {
        return nz_;
    }
    int depth() const {
        return depth_;
    }

    /// Nodes in the lattice, `nx * ny * nz`.
    std::size_t node_count() const {
        return static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) *
               static_cast<std::size_t>(nz_);
    }

private:
    std::size_t index(int i, int j, int k, int q) const {
        return ((static_cast<std::size_t>(i) * static_cast<std::size_t>(ny_) +
                 static_cast<std::size_t>(j)) *
                    static_cast<std::size_t>(nz_) +
                static_cast<std::size_t>(k)) *
                   static_cast<std::size_t>(depth_) +
               static_cast<std::size_t>(q);
    }

    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    int depth_ = 0;
    std::vector<double> data_;
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_FIELD3D_H
