message(STATUS "")
message(STATUS "Looking for dependencies (OpenMP)")

if(WITH_OpenMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP requested but not found: the *_omp programs will be serial.")
    endif()
endif()

add_library(cglbm_dependencies INTERFACE)

if(WITH_OpenMP AND OpenMP_CXX_FOUND)
    target_link_libraries(cglbm_dependencies INTERFACE OpenMP::OpenMP_CXX)
endif()
