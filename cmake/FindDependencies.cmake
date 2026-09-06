message(STATUS "")
message(STATUS "Looking for dependencies (OpenMP, MPI, CUDA)")

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

set(CUDA_FOUND FALSE)
if(WITH_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        set(CMAKE_CUDA_EXTENSIONS OFF)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "CUDA architectures to target")
        endif()
        find_package(CUDAToolkit QUIET)
        set(CUDA_FOUND TRUE)
        message(STATUS "    CUDA ${CMAKE_CUDA_COMPILER_VERSION} found (${CMAKE_CUDA_COMPILER})")
    else()
        message(WARNING "CUDA requested but not found: src/cuda falls back to disabled/host.")
    endif()
endif()

add_library(cglbm_dependencies INTERFACE)

if(WITH_OpenMP AND OpenMP_CXX_FOUND)
    target_link_libraries(cglbm_dependencies INTERFACE OpenMP::OpenMP_CXX)
endif()

# CGLBM_WITH_MPI guards every use of <mpi.h>: without it src/mpi still compiles and behaves as a
# single rank, so the same sources build either way.
if(WITH_MPI AND MPI_CXX_FOUND)
    target_link_libraries(cglbm_dependencies INTERFACE MPI::MPI_CXX)
    target_compile_definitions(cglbm_dependencies INTERFACE CGLBM_WITH_MPI)
    set(CGLBM_MPI_LAUNCHER
        "${MPIEXEC_EXECUTABLE}"
        CACHE FILEPATH "mpirun used by the tests"
    )
endif()

# CGLBM_WITH_CUDA guards every use of CUDA runtime headers and API calls: without it
# src/cuda still compiles and behaves as a disabled host fallback.
if(WITH_CUDA AND CUDA_FOUND)
    target_compile_definitions(cglbm_dependencies INTERFACE CGLBM_WITH_CUDA)
    if(TARGET CUDA::cudart)
        target_link_libraries(cglbm_dependencies INTERFACE CUDA::cudart)
    endif()
endif()
