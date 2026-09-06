message(STATUS "")
message(STATUS "=================================================================================")
message(STATUS "Selecting compiler flags")
message(
    STATUS ${CMAKE_CXX_COMPILER_ID}
           " compiler: "
           ${CMAKE_CXX_COMPILER}
           " (version: "
           ${CMAKE_CXX_COMPILER_VERSION}
           ")"
)
message(STATUS "")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(COMMON_FLAGS
        "-g"
        "-Wall"
        "-Wextra"
        "-Wuninitialized"
        "-Wshadow"
        "-pedantic"
    )
    # The recoloring step divides by the norm of the colour gradient, which goes to zero away from
    # the interface: keep IEEE semantics, never -ffast-math.
    set(RELEASE_FLAGS "-O3" "-fno-fast-math" "-fno-finite-math-only")
    set(DEBUG_FLAGS
        "-Og"
        "-g3"
        "-fsignaling-nans"
        "-ftrapv"
        "-fno-fast-math"
        "-D_GLIBCXX_ASSERTIONS"
    )
    if(ARCH)
        set(COMMON_FLAGS ${COMMON_FLAGS} "-march=${ARCH}")
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(COMMON_FLAGS
        "-g"
        "-Wall"
        "-Wextra"
        "-Wuninitialized"
        "-Wshadow"
        "-pedantic"
    )
    set(RELEASE_FLAGS "-O3" "-fno-fast-math")
    set(DEBUG_FLAGS "-O0" "-g3" "-fno-fast-math" "-D_GLIBCXX_ASSERTIONS")
    if(ARCH)
        set(COMMON_FLAGS ${COMMON_FLAGS} "-march=${ARCH}")
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # icpx compiler (LINUX)
    set(COMMON_FLAGS "-g" "-traceback" "-Wall")
    set(RELEASE_FLAGS "-O3" "-fp-model=precise")
    set(DEBUG_FLAGS "-O0" "-check=stack,uninit" "-fp-model=precise")
    if(ARCH)
        set(COMMON_FLAGS ${COMMON_FLAGS} "-x${ARCH}")
    endif()

else()
    message(WARNING "Unknown compiler ${CMAKE_CXX_COMPILER_ID}: using default flags.")
    set(COMMON_FLAGS "")
    set(RELEASE_FLAGS "-O2")
    set(DEBUG_FLAGS "-O0")
endif()

string(REPLACE ";" " " COMMON_FLAGS_STR "${COMMON_FLAGS}")
string(REPLACE ";" " " RELEASE_FLAGS_STR "${RELEASE_FLAGS}")
string(REPLACE ";" " " DEBUG_FLAGS_STR "${DEBUG_FLAGS}")

if(RELEASE)
    message(STATUS "Release Flags: ${CMAKE_CXX_COMPILER} ${COMMON_FLAGS_STR} ${RELEASE_FLAGS_STR}")
    target_compile_options(cglbm_options_opt INTERFACE ${COMMON_FLAGS} ${RELEASE_FLAGS})
    target_compile_definitions(cglbm_options_opt INTERFACE NDEBUG)
    target_link_options(cglbm_options_opt INTERFACE ${COMMON_FLAGS} ${RELEASE_FLAGS})
endif()
if(DEBUG)
    message(STATUS "Debug Flags: ${CMAKE_CXX_COMPILER} ${COMMON_FLAGS_STR} ${DEBUG_FLAGS_STR}")
    target_compile_options(cglbm_options_dbg INTERFACE ${COMMON_FLAGS} ${DEBUG_FLAGS})
    target_compile_definitions(cglbm_options_dbg INTERFACE DEBUG_CGLBM)
    target_link_options(cglbm_options_dbg INTERFACE ${COMMON_FLAGS} ${DEBUG_FLAGS})
endif()
