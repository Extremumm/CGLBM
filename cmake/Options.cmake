set(_bin_dir_default "${PROJECT_SOURCE_DIR}/bin")
set(BIN_DIR
    "${_bin_dir_default}"
    CACHE PATH "Where to put the binaries"
)
set(_artifacts_dir_default "${PROJECT_SOURCE_DIR}/artifacts")
set(ARTIFACTS_DIR
    "${_artifacts_dir_default}"
    CACHE PATH "Where the tests write their run directories"
)
set(_arch_default "$ENV{ARCH}")
set(ARCH
    "${_arch_default}"
    CACHE STRING "if not empty, value is given to compiler (-march=VALUE or -xVALUE)"
)

set(_release_default ON)
set(_release ${_release_default})
if(DEFINED ENV{RELEASE})
    set(_release "$ENV{RELEASE}")
endif()
option(RELEASE "build optimized version" ${_release})

set(_debug_default ON)
set(_debug ${_debug_default})
if(DEFINED ENV{DEBUG})
    set(_debug "$ENV{DEBUG}")
endif()
option(DEBUG "build debug version" ${_debug})

set(_with_openmp_default ON)
set(_with_openmp ${_with_openmp_default})
if(DEFINED ENV{WITH_OpenMP})
    set(_with_openmp "$ENV{WITH_OpenMP}")
endif()
option(WITH_OpenMP "activate OpenMP" ${_with_openmp})

set(_with_mpi_default ON)
set(_with_mpi ${_with_mpi_default})
if(DEFINED ENV{WITH_MPI})
    set(_with_mpi "$ENV{WITH_MPI}")
endif()
option(WITH_MPI "activate MPI" ${_with_mpi})

set(_with_cuda_default ON)
set(_with_cuda ${_with_cuda_default})
if(DEFINED ENV{WITH_CUDA})
    set(_with_cuda "$ENV{WITH_CUDA}")
endif()
option(WITH_CUDA "activate CUDA" ${_with_cuda})

set(_with_ipo_default OFF)
set(_with_ipo ${_with_ipo_default})
if(DEFINED ENV{WITH_IPO})
    set(_with_ipo "$ENV{WITH_IPO}")
endif()
option(WITH_IPO "Optimization that might be very slow" ${_with_ipo})

set(_with_python_default ON)
set(_with_python ${_with_python_default})
if(DEFINED ENV{WITH_Python})
    set(_with_python "$ENV{WITH_Python}")
endif()
option(WITH_Python "Register the pytest suite with CTest" ${_with_python})

set(_except_default ON)
set(_except ${_except_default})
if(DEFINED ENV{EXCEPT})
    set(_except "$ENV{EXCEPT}")
endif()
option(EXCEPT "Allows you to define per-file exceptions" ${_except})

message(STATUS "")
message(STATUS "   Option        | Default |  Final  | Meaning")
message(STATUS "   ------------- | ------- | ------- | -------")
# Print one formatted row in the option summary table.
function(
    print_cglbm_option
    name
    default
    final
    meaning
)
    set(customized " ")
    if(NOT "${default}" STREQUAL "${final}")
        set(customized "*")
    endif()
    if("${default}" STREQUAL "")
        set(default "''")
    endif()
    if("${final}" STREQUAL "")
        set(final "''")
    endif()
    string(LENGTH "${name}" name_length)
    while(name_length LESS 13)
        string(APPEND name " ")
        math(EXPR name_length "${name_length} + 1")
    endwhile()
    string(LENGTH "${default}" default_length)
    while(default_length LESS 7)
        string(APPEND default " ")
        math(EXPR default_length "${default_length} + 1")
    endwhile()
    string(LENGTH "${final}" final_length)
    while(final_length LESS 7)
        string(APPEND final " ")
        math(EXPR final_length "${final_length} + 1")
    endwhile()
    message(STATUS " ${customized} ${name} | ${default} | ${final} | ${meaning}")
endfunction()

set(cglbm_options
    "RELEASE|${_release_default}|${RELEASE}|Build the optimized version."
    "DEBUG|${_debug_default}|${DEBUG}|Build the debug version."
    "ARCH|${_arch_default}|${ARCH}|Processor architecture. If empty, the result is less optimized."
    "WITH_OpenMP|${_with_openmp_default}|${WITH_OpenMP}|Activate OpenMP."
    "WITH_MPI|${_with_mpi_default}|${WITH_MPI}|Activate MPI."
    "WITH_CUDA|${_with_cuda_default}|${WITH_CUDA}|Activate CUDA."
    "WITH_IPO|${_with_ipo_default}|${WITH_IPO}|Enable link-time/interprocedural optimization."
    "WITH_Python|${_with_python_default}|${WITH_Python}|Register the pytest suite with CTest."
    "EXCEPT|${_except_default}|${EXCEPT}|Apply per-file exceptions."
    "BIN_DIR|${_bin_dir_default}|${BIN_DIR}|Where to put the binaries."
    "ARTIFACTS_DIR|${_artifacts_dir_default}|${ARTIFACTS_DIR}|Where the tests write their runs."
)
foreach(option_entry IN LISTS cglbm_options)
    string(REPLACE "|" ";" option_fields "${option_entry}")
    list(GET option_fields 0 option_name)
    list(GET option_fields 1 option_default)
    list(GET option_fields 2 option_final)
    list(GET option_fields 3 option_meaning)
    print_cglbm_option("${option_name}" "${option_default}" "${option_final}" "${option_meaning}")
endforeach()
