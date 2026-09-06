#include "lbm/output_writer.h"

#include <stdexcept>
#include <string>

namespace cglbm {
namespace lbm {
namespace {

/// Close `file` and complain if anything went wrong along the way.
///
/// A stream only reports a write failure at the next check, and buffered data
/// is not flushed until close, so both are looked at here rather than after
/// each `<<`.
void close_checked(std::ofstream& file, const std::string& filename) {
    file.close();
    if (!file) {
        throw OutputError("failed to write " + filename);
    }
}

}  // namespace

std::ofstream CsvWriter::open_checked(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw OutputError("cannot open " + filename + " for writing");
    }
    file.precision(precision_);
    return file;
}

void CsvWriter::write_grids(int timestep, const MacroscopicState& state) const {
    const std::string suffix = "_" + std::to_string(timestep) + ".csv";
    const std::string density_name = "density" + suffix;
    const std::string velocity_name = "velocity" + suffix;
    const std::string phase_name = "phase" + suffix;
    const std::string pressure_name = "pressure" + suffix;

    std::ofstream density = open_checked(density_name);
    std::ofstream velocity = open_checked(velocity_name);
    std::ofstream phase = open_checked(phase_name);
    std::ofstream pressure = open_checked(pressure_name);

    const int nx = state.density->nx();
    const int ny = state.density->ny();

    // Row j of the file is lattice row j, so the grids load as array[y, x].
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            density << (*state.density)(i, j);
            velocity << (*state.velocity)(i, j, 0) << "," << (*state.velocity)(i, j, 1);
            phase << (*state.phase)(i, j);
            pressure << (*state.pressure)(i, j);
            if (i < nx - 1) {
                density << ",";
                velocity << ",";
                phase << ",";
                pressure << ",";
            }
        }
        density << "\n";
        velocity << "\n";
        phase << "\n";
        pressure << "\n";
    }

    close_checked(density, density_name);
    close_checked(velocity, velocity_name);
    close_checked(phase, phase_name);
    close_checked(pressure, pressure_name);
}

void CsvWriter::open_interface_track() {
    interface_file_.open("interface.csv");
    if (!interface_file_) {
        throw OutputError("cannot open interface.csv for writing");
    }
    interface_file_.precision(precision_);
    interface_file_ << "Timestep, phi, y" << std::endl;
}

void CsvWriter::write_interface(int timestep, const Field& phase) {
    const int column = phase.nx() / 2;
    for (int j = 0; j < phase.ny(); ++j) {
        const double phi_local = phase(column, j);
        if (phi_local > -1.0 && phi_local < 1.0) {
            interface_file_ << timestep << ", " << phi_local << ", " << j << std::endl;
        }
    }
    if (!interface_file_) {
        throw OutputError("failed to write interface.csv");
    }
}

void write_vtk(const std::string& prefix, int timestep, const MacroscopicState& state) {
    const std::string filename = prefix + std::to_string(timestep) + ".vtk";
    std::ofstream file(filename);
    if (!file) {
        throw OutputError("cannot open " + filename + " for writing");
    }

    const int nx = state.density->nx();
    const int ny = state.density->ny();

    file << "# vtk DataFile Version 4.0\n"
         << "Lattice Boltzmann Method Data\n"
         << "ASCII\n"
         << "DATASET STRUCTURED_POINTS\n"
         << "DIMENSIONS " << nx << " " << ny << " 1\n"
         << "ORIGIN 0 0 0\n"
         << "SPACING 1 1 1\n";

    file << "POINT_DATA " << nx * ny << "\n";
    file << "SCALARS density float\nLOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << (*state.density)(i, j) << "\n";
        }
    }

    file << "VECTORS velocity float\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << (*state.velocity)(i, j, 0) << " " << (*state.velocity)(i, j, 1) << " 0.0\n";
        }
    }

    file << "SCALARS phase_field float\nLOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << (*state.phase)(i, j) << "\n";
        }
    }

    file << "SCALARS pressure float\nLOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << (*state.pressure)(i, j) << "\n";
        }
    }

    close_checked(file, filename);
}

}  // namespace lbm
}  // namespace cglbm
