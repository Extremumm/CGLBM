message(STATUS "")
message(STATUS "=================================================================================")
message(STATUS "Preparing compilation of CGLBM")
message(STATUS "")
message(STATUS "If you need to change compilation option, change:")
message(STATUS " - the cmake/CompilerFlags.cmake for every file")
message(STATUS " - the cmake/Exceptions.cmake for each file (if necessary)")
message(STATUS "")
message(STATUS "If you need to adapt for a new computer please:")
message(STATUS " - adapt the dependencies in cmake/FindDependencies.cmake")
message(STATUS " - fix the compiler flags in cmake/CompilerFlags.cmake")
message(STATUS "")
message(STATUS "=================================================================================")

cmake_minimum_required(VERSION 3.16)

project(cglbm CXX)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds not allowed. \
            Please make a new directory (called a build directory) and run CMake from there.\n"
    )
endif()

# OPTIONS ###
include(Options)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(WITH_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(
        RESULT ipo_supported
        OUTPUT ipo_error
        LANGUAGES CXX
    )
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "IPO not supported")
    endif()
endif()

# Dependencies ###
include(FindDependencies)

# definitions of the targets ###
set(opt_versions "")
if(RELEASE)
    list(APPEND opt_versions "opt")
    add_custom_target(opt COMMENT "Build all optimized targets")
    add_library(cglbm_options_opt INTERFACE)
endif()
if(DEBUG)
    list(APPEND opt_versions "dbg")
    add_custom_target(dbg COMMENT "Build all debug targets")
    add_library(cglbm_options_dbg INTERFACE)
endif()
if(NOT opt_versions)
    message(FATAL_ERROR "Nothing to build: enable at least one of RELEASE and DEBUG.")
endif()

file(
    GLOB_RECURSE
    program_paths
    CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/programs/*/main_*.cpp"
)
# append a fictive sourcefile for producing the library
set(program_paths "${PROJECT_SOURCE_DIR}/src/main_cglbm.cpp" ${program_paths})
list(REMOVE_DUPLICATES program_paths)

message(STATUS "")
message(STATUS "=================================================================================")
foreach(program_path IN LISTS program_paths)
    get_filename_component(program_dir ${program_path} DIRECTORY)
    get_filename_component(base_target_name ${program_path} NAME_WE)
    string(REGEX REPLACE "^main_" "" base_target_name ${base_target_name})

    # every non-main source living next to the main_*.cpp belongs to that program
    file(
        GLOB_RECURSE
        SRC_FILES
        CONFIGURE_DEPENDS
        ${program_dir}/[!_]*.cpp
    )

    add_custom_target(
        "${base_target_name}_all"
        COMMENT "Build both ${base_target_name}_opt and ${base_target_name}_dbg target"
    )

    if(base_target_name STREQUAL "cglbm")
        # The library is the whole of src/ minus the fictive entry point. It is
        # still header-only: the solver lives in the monolithic programs until
        # src/lbm gains its implementation units.
        list(
            FILTER
            SRC_FILES
            EXCLUDE
            REGEX
            "main_cglbm\\.cpp$"
        )
    endif()

    foreach(opt_version IN LISTS opt_versions)
        set(target_name "${base_target_name}_${opt_version}")

        set(target_is_interface FALSE)
        if(base_target_name STREQUAL "cglbm")
            if(SRC_FILES)
                add_library(${target_name} ${SRC_FILES})
                target_include_directories(${target_name} PUBLIC "${PROJECT_SOURCE_DIR}/src")
                target_link_libraries(
                    ${target_name} PRIVATE cglbm_dependencies cglbm_options_${opt_version}
                )
            else()
                # header-only for now: an INTERFACE target carries no build rule
                set(target_is_interface TRUE)
                add_library(${target_name} INTERFACE)
                target_include_directories(${target_name} INTERFACE "${PROJECT_SOURCE_DIR}/src")
            endif()
            add_library(cglbm_program_${opt_version} INTERFACE)
            target_link_libraries(cglbm_program_${opt_version} INTERFACE ${target_name})
        else()
            add_executable(${target_name} ${SRC_FILES})
            target_link_libraries(
                ${target_name} PRIVATE cglbm_dependencies cglbm_options_${opt_version}
                                       cglbm_program_${opt_version}
            )
            # mirror the programs/ tree inside bin/, as documented in the README
            file(RELATIVE_PATH program_subdir "${PROJECT_SOURCE_DIR}/programs" "${program_dir}")
            set_target_properties(
                ${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}/${program_subdir}"
            )
        endif()

        if(NOT target_is_interface)
            add_dependencies("${base_target_name}_all" ${target_name})
            add_dependencies("${opt_version}" ${target_name})
        endif()
    endforeach(opt_version)

    if(NOT base_target_name STREQUAL "cglbm")
        message(STATUS " - ${base_target_name}")
    endif()
endforeach(program_path)

include(CompilerFlags)
include(Exceptions)

message(STATUS "")
message(STATUS "=================================================================================")

# Computer specific ###
include(Tests)
