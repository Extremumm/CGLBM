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

const char* interface_viscosity_name(InterfaceViscosity mixing) {
    switch (mixing) {
    case InterfaceViscosity::Harmonic:
        return "harmonic";
    case InterfaceViscosity::Laminate:
        return "laminate";
    case InterfaceViscosity::Arithmetic:
    default:
        return "arithmetic";
    }
}

bool take_viscosity_option(int* argc, char** argv, InterfaceViscosity* mixing) {
    const std::string prefix = "--viscosity=";
    bool valid = true;
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        const std::string argument = argv[n];
        if (argument.rfind(prefix, 0) != 0) {
            argv[kept++] = argv[n];
            continue;
        }
        const std::string value = argument.substr(prefix.size());
        if (value == "arithmetic") {
            *mixing = InterfaceViscosity::Arithmetic;
        } else if (value == "harmonic") {
            *mixing = InterfaceViscosity::Harmonic;
        } else if (value == "laminate") {
            *mixing = InterfaceViscosity::Laminate;
        } else {
            valid = false;
        }
    }
    *argc = kept;
    return valid;
}

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm
