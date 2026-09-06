if(WITH_Python)
    include(CTest)

    find_package(Python3 3.9 COMPONENTS Interpreter)
    if(NOT Python3_FOUND)
        message(WARNING "Python3 not found: Python tests will not be configured.")
        return()
    endif()

    # Prefer the pytest of the interpreter we found; fall back to a pytest on PATH, so that a
    # virtualenv and a system install both work.
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import pytest"
        RESULT_VARIABLE Pytest_IMPORT_RESULT
        OUTPUT_QUIET ERROR_QUIET
    )
    if(Pytest_IMPORT_RESULT EQUAL 0)
        set(PYTEST_COMMAND ${Python3_EXECUTABLE} -m pytest)
    else()
        find_program(Pytest_EXECUTABLE pytest)
        if(NOT Pytest_EXECUTABLE)
            message(WARNING "pytest not found for ${Python3_EXECUTABLE}: "
                            "Python tests will not be configured (pip install -r requirements.txt)."
            )
            return()
        endif()
        set(PYTEST_COMMAND ${Pytest_EXECUTABLE})
    endif()

    # pytest exits with 5 when a file collects no test for the selected marker, which is not a
    # failure for us: every file is run once per marker.
    string(
        CONCAT PYTEST_RUNNER
               "import subprocess\n"
               "import sys\n"
               "code = subprocess.run(sys.argv[1:]).returncode\n"
               "sys.exit(0 if code == 5 else code)\n"
    )

    file(
        GLOB_RECURSE CGLBM_PYTEST_FILES CONFIGURE_DEPENDS
        RELATIVE "${PROJECT_SOURCE_DIR}/programs"
        "${PROJECT_SOURCE_DIR}/programs/test_*.py"
    )

    message(STATUS "")
    message(STATUS "Registering the pytest suite:")
    foreach(pytest_file IN LISTS CGLBM_PYTEST_FILES)
        string(REPLACE "/" "." pytest_stem "${pytest_file}")
        string(REGEX REPLACE "\\.py$" "" pytest_stem "${pytest_stem}")

        foreach(pytest_marker unit_test validation verification)
            add_test(
                NAME "cglbm-${pytest_marker}-${pytest_stem}"
                COMMAND ${Python3_EXECUTABLE} -c "${PYTEST_RUNNER}" ${PYTEST_COMMAND}
                        "${PROJECT_SOURCE_DIR}/programs/${pytest_file}" -m "${pytest_marker}" -v
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            )
            set_tests_properties(
                "cglbm-${pytest_marker}-${pytest_stem}"
                PROPERTIES ENVIRONMENT
                           "CGLBM_BIN_DIR=${BIN_DIR};CGLBM_ARTIFACTS_DIR=${ARTIFACTS_DIR}"
            )
        endforeach()
        message(STATUS " - ${pytest_file}")
    endforeach()
endif()
