/// Naming device for the library target, not a program.
///
/// CMakeLists.txt discovers programs by globbing `main_*.cpp` and takes the
/// target name from the file. This file is what makes the whole of `src/`
/// build as `libcglbm_opt` / `libcglbm_dbg`; CMake filters it back out of the
/// library's sources, so nothing here is ever compiled.
///
/// The solvers are the programs under `programs/solvers`. Each is a case
/// definition handed to `cglbm::lbm::Solver`, which is where the scheme lives.

int main() {
    return 0;
}
