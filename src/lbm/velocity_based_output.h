#ifndef CGLBM_LBM_VELOCITY_BASED_OUTPUT_H
#define CGLBM_LBM_VELOCITY_BASED_OUTPUT_H

#include "lbm/velocity_based_solver.h"

/// What the velocity-based programs share: the fields they write, and how
/// they read a number from their command line.
namespace cglbm {
namespace lbm {
namespace velocity_based {

/// density_<t>.csv, velocity_<t>.csv (u_x,u_y per node), phase_<t>.csv and
/// pressure_<t>.csv in the working directory: a row per j, a column per i,
/// ten significant digits, the layout the colour-gradient programs write and
/// pycglbm.CaseOutput reads. The phase is psi = 2c - 1, the pressure relative
/// to the far field.
void write_fields(const Solver& solver, int timestep);

/// The whole of `text` as a finite number; false, with `value` untouched,
/// when it is not one.
bool parse_number(const char* text, double* value);

/// Name of the mixing, as take_viscosity_option reads it: "arithmetic",
/// "harmonic" or "laminate".
const char* interface_viscosity_name(InterfaceViscosity mixing);

/// Takes every --viscosity=arithmetic|harmonic|laminate out of argv, leaving
/// the other arguments in order, and sets `mixing` from the last valid one.
/// False if any of them names none of the three.
bool take_viscosity_option(int* argc, char** argv, InterfaceViscosity* mixing);

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_VELOCITY_BASED_OUTPUT_H
