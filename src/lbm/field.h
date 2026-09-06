#ifndef CGLBM_LBM_FIELD_H
#define CGLBM_LBM_FIELD_H

#include <cstddef>
#include <vector>

/// Lattice-shaped storage.
///
/// The solvers used to declare their lattice as file-scope arrays,
/// `double f[Lx][Ly][Q]`, which put the whole state in `.bss` and fixed the
/// domain size at compile time. At the production resolution of the
/// Rayleigh-Taylor case that came to about 2 GB of static data and needed
/// `-mcmodel=medium` to link at all -- see `cmake/Exceptions.cmake`.
///
/// `Field` holds the same values in one heap allocation, so the lattice size
/// becomes a runtime parameter and the binary stays small. The memory layout is
/// unchanged: element (i, j, k) sits at `(i * ny + j) * depth + k`, which is
/// what `[Lx][Ly][Q]` gave and what `isotropic_gradient.h` expects of a
/// `depth == 1` field.
///
/// Storage is value-initialised to zero, matching the `.bss` the arrays came
/// from. Two solvers rely on that: the wall-bounded cases never write the
/// source term on the `j == 0` row and read the zeros back.

namespace cglbm {
namespace lbm {

class Field {
  public:
    Field() = default;

    /// `nx` by `ny` nodes, `depth` values each, all zero.
    Field(int nx, int ny, int depth = 1)
        : nx_(nx), ny_(ny), depth_(depth),
          data_(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                    static_cast<std::size_t>(depth),
                0.0) {}

    /// Value at node (i, j) of a scalar field.
    double& operator()(int i, int j) { return data_[index(i, j, 0)]; }
    double operator()(int i, int j) const { return data_[index(i, j, 0)]; }

    /// Component k at node (i, j).
    double& operator()(int i, int j, int k) { return data_[index(i, j, k)]; }
    double operator()(int i, int j, int k) const { return data_[index(i, j, k)]; }

    /// Contiguous storage, laid out as `[nx][ny][depth]`.
    ///
    /// A scalar field's pointer is what the gradient stencils take.
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int depth() const { return depth_; }

  private:
    std::size_t index(int i, int j, int k) const {
        return (static_cast<std::size_t>(i) * static_cast<std::size_t>(ny_) +
                static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(depth_) +
               static_cast<std::size_t>(k);
    }

    int nx_ = 0;
    int ny_ = 0;
    int depth_ = 0;
    std::vector<double> data_;
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_FIELD_H
