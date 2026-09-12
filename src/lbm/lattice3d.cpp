#include "lbm/lattice3d.h"

#include <cstring>

namespace cglbm {
namespace lbm {

bool lattice_3d_from_name(const char* name, Lattice3DKind* kind) {
    if (std::strcmp(name, "D3Q19") == 0 || std::strcmp(name, "d3q19") == 0) {
        *kind = Lattice3DKind::D3Q19;
        return true;
    }
    if (std::strcmp(name, "D3Q27") == 0 || std::strcmp(name, "d3q27") == 0) {
        *kind = Lattice3DKind::D3Q27;
        return true;
    }
    return false;
}

const char* lattice_3d_name(Lattice3DKind kind) {
    return kind == Lattice3DKind::D3Q27 ? "D3Q27" : "D3Q19";
}

}  // namespace lbm
}  // namespace cglbm
