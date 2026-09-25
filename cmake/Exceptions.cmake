if(EXCEPT)
    message(STATUS "")
    message(STATUS "=============================================================================")
    message(STATUS "Applying some exceptions:")

    # No per-file exception is needed at the moment.
    #
    # The OpenMP Rayleigh-Taylor program used to need -mcmodel=medium: it declared its lattice as
    # static arrays of 1024 x 4096 x 9 doubles, about 2 GB of .bss, which overflows the 32-bit
    # displacements of the default small code model on x86-64 ELF and fails at link time with
    # "relocation truncated to fit". src/lbm/field.h puts the lattice on the heap, so the binary is
    # small again and the exception is gone.
    message(STATUS " - none")

endif()
