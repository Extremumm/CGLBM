#include "lbm/velocity_based_output.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

namespace cglbm {
namespace lbm {
namespace velocity_based {

void write_fields(const Solver& solver, int timestep) {
    const std::string suffix = "_" + std::to_string(timestep) + ".csv";
    std::ofstream density("density" + suffix);
    std::ofstream velocity("velocity" + suffix);
    std::ofstream phase("phase" + suffix);
    std::ofstream pressure("pressure" + suffix);
    for (auto* file : {&density, &velocity, &phase, &pressure}) {
        file->precision(10);
    }
    const int nx = solver.nx();
    for (int j = 0; j < solver.ny(); ++j) {
        for (int i = 0; i < nx; ++i) {
            const char* separator = i < nx - 1 ? "," : "\n";
            density << solver.density(i, j) << separator;
            velocity << solver.velocity_x(i, j) << "," << solver.velocity_y(i, j) << separator;
            phase << solver.phase(i, j) << separator;
            pressure << solver.pressure(i, j) << separator;
        }
    }
}

bool parse_number(const char* text, double* value) {
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm
