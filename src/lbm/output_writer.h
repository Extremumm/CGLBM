#ifndef CGLBM_LBM_OUTPUT_WRITER_H
#define CGLBM_LBM_OUTPUT_WRITER_H

#include <fstream>
#include <stdexcept>
#include <string>

#include "lbm/field.h"

/// Writing a run to disk.
///
/// Every program writes into the current working directory under fixed names,
/// so a run needs a directory of its own. That contract is unchanged; what is
/// new is that a failure to write is reported. The solvers used to call
/// `ofstream::open` and never look at the result, so a full disk or a
/// read-only directory produced a silent, empty run.

namespace cglbm {
namespace lbm {

/// The macroscopic state at one instant, as the writers see it.
struct MacroscopicState {
    const Field* density;   ///< mixture density, one value per node
    const Field* velocity;  ///< velocity, two components per node
    const Field* phase;     ///< phase field in [-1, 1]
    const Field* pressure;  ///< mixture pressure
};

/// Thrown when an output file cannot be opened or written.
///
/// A run that cannot record its result has failed, and a solver that keeps
/// going wastes the hours it still has to run.
class OutputError : public std::runtime_error {
public:
    explicit OutputError(const std::string& message) : std::runtime_error(message) {}
};

/// Writes the CSV grids and, optionally, the interface track.
///
/// The four grids are `density_<t>.csv`, `velocity_<t>.csv`, `phase_<t>.csv`
/// and `pressure_<t>.csv`: one row per lattice row j, one column per node i,
/// so `numpy` reads them as `array[y, x]`. Velocity interleaves its two
/// components, giving `ny` rows of `2 * nx` values.
class CsvWriter {
public:
    /// `precision` is the number of significant digits; 6 is the iostream
    /// default and what every run before this class wrote.
    explicit CsvWriter(int precision) : precision_(precision) {}

    /// Write the four grids for `timestep`.
    void write_grids(int timestep, const MacroscopicState& state) const;

    /// Start `interface.csv` and write its header.
    void open_interface_track();

    /// Append the interface crossing along the mid-plane, i = nx / 2.
    ///
    /// One line per node whose phase field is strictly inside [-1, 1], which is
    /// the diffuse interface. Called every step, not every interval: the
    /// capillary cases measure an oscillation period from it.
    void write_interface(int timestep, const Field& phase);

private:
    std::ofstream open_checked(const std::string& filename) const;

    int precision_;
    std::ofstream interface_file_;
};

/// Write the whole state as one legacy-ASCII VTK file, for ParaView or VisIt.
///
/// `prefix` is completed with the timestep and `.vtk`.
void write_vtk(const std::string& prefix, int timestep, const MacroscopicState& state);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_OUTPUT_WRITER_H
