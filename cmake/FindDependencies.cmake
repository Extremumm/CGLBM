message(STATUS "")
message(STATUS "Looking for dependencies (OpenMP, MPI)")

if(WITH_OpenMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP requested but not found: the *_omp programs will be serial.")
    endif()
endif()

if(WITH_MPI)
    find_package(MPI QUIET COMPONENTS CXX)
    if(NOT MPI_CXX_FOUND)
        message(WARNING "MPI requested but not found: src/mpi falls back to a single rank.")
    else()
        message(STATUS "    MPI ${MPI_CXX_VERSION} found")
    endif()
endif()

add_library(cglbm_dependencies INTERFACE)

if(WITH_OpenMP AND OpenMP_CXX_FOUND)
    target_link_libraries(cglbm_dependencies INTERFACE OpenMP::OpenMP_CXX)
endif()

# CGLBM_WITH_MPI guards every use of <mpi.h>: without it src/mpi still compiles
# and behaves as a single rank, so the same sources build either way.
if(WITH_MPI AND MPI_CXX_FOUND)
    target_link_libraries(cglbm_dependencies INTERFACE MPI::MPI_CXX)
    target_compile_definitions(cglbm_dependencies INTERFACE CGLBM_WITH_MPI)
    set(CGLBM_MPI_LAUNCHER "${MPIEXEC_EXECUTABLE}" CACHE FILEPATH "mpirun used by the tests")
endif()
