# Storage and lattices

## `Field` — `src/lbm/field.h`

Lattice-shaped storage of `double`, one heap allocation, value-initialised to
zero.

```cpp
class Field {
public:
    Field() = default;                              // empty: nx = ny = depth = 0
    Field(int nx, int ny, int depth = 1);           // nx * ny * depth zeros
    double& operator()(int i, int j);               // component 0 of node (i, j)
    double& operator()(int i, int j, int k);        // component k of node (i, j)
    double* data();                                 // contiguous [nx][ny][depth]
    int nx() const; int ny() const; int depth() const;
};
```

Element `(i, j, k)` is at `(i * ny + j) * depth + k`. Indices are not checked.
`operator()(i, j)` is `operator()(i, j, 0)`, so a scalar field and component 0
of a vector field read the same way. `data()` of a scalar field is the
`[i * ny + j]` array the gradient stencils of
[`isotropic_gradient.h`](numerical-kernels.md#gradient-stencils) take.

How the solvers use `depth`:

| Field | depth |
|---|---|
| scalar (density, pressure, phase, …) | 1 |
| velocity, force | 2 |
| populations `f`, `g`, `f1`, `f2`, operators `omega_*`, source | `kQ` = 9 |

## `Field3D` — `src/lbm/field3d.h`

The same with one more axis. Element `(i, j, k, q)` is at
`((i * ny + j) * nz + k) * depth + q`.

```cpp
class Field3D {
public:
    Field3D(int nx, int ny, int nz, int depth = 1);
    double& operator()(int i, int j, int k);          // component 0
    double& operator()(int i, int j, int k, int q);   // component q
    double* data();
    int nx() const; int ny() const; int nz() const; int depth() const;
    std::size_t node_count() const;                   // nx * ny * nz
};
```

It is a separate class so that the 2D call `f(i, j, k)` with `k` the velocity
index cannot be mistaken for a 3D node.

## D2Q9 — `src/lbm/d2q9.h`

```
    6   2   5
     \  |  /
    3 - 0 - 1
     /  |  \
    7   4   8
```

| `k` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| `kXi[k]` | (0,0) | (1,0) | (0,1) | (−1,0) | (0,−1) | (1,1) | (−1,1) | (−1,−1) | (1,−1) |
| `kW[k]` | 4/9 | 1/9 | 1/9 | 1/9 | 1/9 | 1/36 | 1/36 | 1/36 | 1/36 |

The half-way bounce-back of the colour-gradient solvers relies on this order:
the directions leaving through the bottom wall, 4, 7, 8, reflect into
`k − 2` = 2, 5, 6, and those leaving through the top wall, 2, 5, 6, into
`k + 2`.

Other constants in the header:

| Name | Value | Meaning |
|---|---|---|
| `kQ` | 9 | number of velocities |
| `kGradientEpsilon` | 1e-10 | a colour-gradient norm below this is treated as zero before dividing by it (Ω⁽²⁾, recolouring) |
| `kInterfaceGradientFloor` | 1e-6 | a node whose `|∇φ_N|` is below this carries no interface: its normal, curvature and capillary force are zero |

`LatticeUnits` holds `dx = 1` and `dt = 1` and derives the sound speed:

```cpp
struct LatticeUnits {
    double dx = 1.0, dt = 1.0;
    double cs()  const;   // dx / (sqrt(3) dt)
    double cs2() const;   // cs^2 = 1/3
    double cs4() const;   // cs^4
    double cs6() const;   // cs^6
};
```

The velocity-based scheme has its own copy of the velocity set, in the same
order, in `velocity_based.h` (`kVelocity`, `kWeight`, `kOpposite`,
`kSoundSpeedSquared = 1/3`).

## D3Q19 — `src/lbm/d3q19.h`

`kQ3D = 19`: rest, the six axial directions, the twelve edge directions.

| `q` | 0 | 1–6 | 7–18 |
|---|---|---|---|
| `kXi3D[q]` | (0,0,0) | (±1,0,0), (0,±1,0), (0,0,±1), in the order +x, −x, +y, −y, +z, −z | (1,1,0), (−1,−1,0), (1,−1,0), (−1,1,0), (1,0,1), (−1,0,−1), (1,0,−1), (−1,0,1), (0,1,1), (0,−1,−1), (0,1,−1), (0,−1,1) |
| `kW3D[q]` | 1/3 | 1/18 | 1/36 |

`kOpposite3D[q]` is the index of `−kXi3D[q]`: `{0, 2, 1, 4, 3, 6, 5, 8, 7,
10, 9, 12, 11, 14, 13, 16, 15, 18, 17}`. Opposite directions are adjacent
pairs. `speed_squared_3d(q)` returns `|kXi3D[q]|²`.

## D3Q27 — `src/lbm/d3q27.h`

`kQ3D27 = 27`. The first nineteen velocities are those of D3Q19 in the same
order; entries 19–26 are the eight corners `(±1, ±1, ±1)`, also in opposite
pairs. Weights `kW3D27`: 8/27 at rest, 2/27 axial, 1/54 edge, 1/216 corner.
`kOpposite3D27` maps each velocity to its reverse. The tables are written one
row per line and fenced with `// clang-format off` so that the layout survives
the formatter.

## `Lattice3D` — `src/lbm/lattice3d.h`

One descriptor per 3D lattice, so that `TwoPopulationSolver3D` has a single
body for both.

```cpp
enum class Lattice3DKind { D3Q19, D3Q27 };

struct Lattice3D {
    const char* name;               // "D3Q19" or "D3Q27"
    int q;                          // 19 or 27
    const double (*xi)[3];          // velocities
    const double* w;                // weights of the standard equilibrium
    const int* opposite;            // reversed direction, for bounce-back
    double cs2_coefficient;         // (c_s^k)^2 = cs2_coefficient * (1 - alpha_k)
    double shell_coefficient[4];    // phi_q^k = shell_coefficient[|e_q|^2] * (1 - alpha_k), q != 0
};
```

| | `cs2_coefficient` | `shell_coefficient` (indexed by `|e|²` = 0, 1, 2, 3) |
|---|---|---|
| `kLatticeD3Q19` | 1/2 | {0, 1/12, 1/24, 0} |
| `kLatticeD3Q27` | 9/19 | {0, 2/19, 1/38, 1/152} |

These are the solutions of the rest-weight conditions of the two-population
model: `Σ_q φ_q = 1`, fourth-order isotropy, and on D3Q27 one sixth-order
relation. Entry 0 is unused because the rest weight is `α_k` itself.

Functions:

| Function | Returns |
|---|---|
| `lattice_3d(kind)` | `kLatticeD3Q27` for `D3Q27`, `kLatticeD3Q19` otherwise |
| `speed_squared_3d(lattice, q)` | `|e_q|²`, 0 to 3 |
| `rest_weight_3d(lattice, q, alpha)` | `alpha` for `q = 0`, else `shell_coefficient[|e_q|²] · (1 − alpha)` |
| `sound_speed_squared_3d(lattice, alpha)` | `cs2_coefficient · (1 − alpha)` |
| `lattice_3d_from_name(name, &kind)` | true and sets `kind` for exactly `"D3Q19"`, `"d3q19"`, `"D3Q27"` or `"d3q27"`; false otherwise, `kind` untouched |
| `lattice_3d_name(kind)` | `"D3Q27"` or `"D3Q19"` |

The amplitude of the enhanced equilibrium is deliberately *not* stored here;
`TwoPopulationSolver3D` derives it from the weights (see
[Two-population solvers](two-population-solvers.md#the-3d-equilibrium)).
