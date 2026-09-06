if(EXCEPT)
    message(STATUS "")
    message(STATUS "=============================================================================")
    message(STATUS "Applying some exceptions:")

    # The OpenMP Rayleigh-Taylor program declares its lattice as static arrays of 1024 x 4096 x 9
    # doubles, i.e. about 2 GB of .bss. That overflows the 32-bit displacements of the default small
    # code model on x86-64 ELF and fails at link time with "relocation truncated to fit". The medium
    # code model puts the large objects in a far data section.
    foreach(opt_version IN LISTS opt_versions)
        set(_target "rayleigh_taylor_omp_${opt_version}")
        if(TARGET ${_target}
           AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64"
           AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"
        )
            message(STATUS " - ${_target}: -mcmodel=medium (~2 GB of static lattice arrays)")
            target_compile_options(${_target} PRIVATE "-mcmodel=medium")
            target_link_options(${_target} PRIVATE "-mcmodel=medium")
        endif()
    endforeach()

endif()
