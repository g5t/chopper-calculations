cmake_minimum_required(VERSION 3.26)
project(chopper-calculations)

set(CMAKE_CXX_STANDARD 20)

# Warn if the user invokes CMake directly
if (NOT SKBUILD)
    message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

find_package(Python 3.8 COMPONENTS
        REQUIRED COMPONENTS Interpreter Development.Module
        OPTIONAL_COMPONENTS Development.SABIModule
        )

if (NOT SKBUILD)
    execute_process(
            COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
endif()

# An installed chopper-lib -- conda-forge's, say -- if there is one, otherwise the tagged
# source from GitHub, built static and linked into the extensions so a wheel carries no
# library of its own.
include(FetchContent)
set(BUILD_SHARED_LIBS OFF)
FetchContent_Declare(chopper_lib
        GIT_REPOSITORY https://github.com/mcdotstar/mcstas-chopper-lib.git
        GIT_TAG v4.4.2
        FIND_PACKAGE_ARGS 4.4.2 CONFIG)
FetchContent_MakeAvailable(chopper_lib)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
        # Name of the extension
        _chopcal_impl

        # Target the stable ABI for Python 3.12+, which reduces
        # the number of binary wheels that must be built. This
        # does nothing on older Python versions
        STABLE_ABI

        # Source code goes here
        src/chopcal.cpp
        src/choppers.cpp
)
target_link_libraries(_chopcal_impl PRIVATE chopper_lib::chopper_lib)

nanobind_add_module(_chopper_lib_impl
        STABLE_ABI
        src/mcstas_chopper_lib.cpp
)
target_link_libraries(_chopper_lib_impl PRIVATE chopper_lib::chopper_lib)

# The bindings convert between wavelength and inverse velocity, and have to do it with
# the same numbers the library does or the two disagree in the eighth digit. The library
# was compiled with CHOPPER_LIB_DEFINITIONS -- the McStas runtime's V2K, K2V and PI, set
# by its package file or left in the cache by FetchContent -- so re-emit them under
# prefixed names the C++ can read without colliding with the chopcal::constants
# variables named PI, V2K and K2V. Those are derived from SI/CODATA; agreeing with McStas
# is the point of these.
set(CHOPCAL_RUNTIME_DEFINITIONS)
foreach (definition IN LISTS CHOPPER_LIB_DEFINITIONS)
    if (definition MATCHES "^([A-Za-z0-9_]+)=(.*)$")
        list(APPEND CHOPCAL_RUNTIME_DEFINITIONS
                "CHOPPER_LIB_${CMAKE_MATCH_1}=${CMAKE_MATCH_2}")
    endif ()
endforeach ()
target_compile_definitions(_chopper_lib_impl PRIVATE ${CHOPCAL_RUNTIME_DEFINITIONS})

# Install directive for scikit-build-core
install(TARGETS _chopcal_impl LIBRARY DESTINATION chopcal)
install(TARGETS _chopper_lib_impl LIBRARY DESTINATION chopcal)

#if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
#    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#endif()
#
## Detect the installed nanobind package and import it into CMake
#execute_process(
#        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
#        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
#list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
#find_package(nanobind CONFIG REQUIRED)
#
#add_subdirectory(src)
