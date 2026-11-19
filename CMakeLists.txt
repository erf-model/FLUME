# ---------------------------------------------------------------------------
# FLUME -- top-level CMake build (stand-alone)
#
#   git submodule update --init Submodules/ERF
#   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX
#   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-rel -j 8
#   ctest --test-dir build-rel --output-on-failure
#
# AMReX is the one ERF pins: Submodules/ERF/Submodules/AMReX. FLUME has no
# AMReX pin of its own. The stand-alone build compiles no ERF code.
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.25)   # required by the pinned AMReX

project(flume
        VERSION 0.1.0
        DESCRIPTION "Two-dimensional shallow-water flood solver on AMReX"
        LANGUAGES C CXX)

# AMReX at the pinned commit uses C++20 in its headers.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(FLUME_MPI          "Build with MPI"                                       ON)
option(FLUME_OMP          "Build with OpenMP"                                    OFF)
option(FLUME_TESTING      "AMReX_TESTING + assertions (SNaN fills, as GPU CI)"   OFF)
option(FLUME_ENABLE_TESTS "Register the tests"                                   ON)
option(FLUME_WERROR       "Warnings in FLUME's own sources are errors"           OFF)

set(FLUME_PRECISION DOUBLE CACHE STRING "Floating-point precision: DOUBLE or SINGLE")
set_property(CACHE FLUME_PRECISION PROPERTY STRINGS DOUBLE SINGLE)

# GPU backend, passed to AMReX. HIP and SYCL need the vendor compiler as
# CMAKE_CXX_COMPILER (ROCm clang++, icpx) and, for HIP, AMReX_AMD_ARCH.
set(FLUME_GPU_BACKEND NONE CACHE STRING "GPU backend: NONE, CUDA, HIP or SYCL")
set_property(CACHE FLUME_GPU_BACKEND PROPERTY STRINGS NONE CUDA HIP SYCL)

# Extra flags for mpiexec in the tests (CI passes --oversubscribe to Open MPI
# so 4-rank tests run on runners with fewer cores).
set(FLUME_MPIEXEC_PREFLAGS "" CACHE STRING "Extra mpiexec flags for the tests")

foreach(opt PRECISION GPU_BACKEND)
    get_property(allowed CACHE FLUME_${opt} PROPERTY STRINGS)
    if(NOT FLUME_${opt} IN_LIST allowed)
        message(FATAL_ERROR "FLUME_${opt}=${FLUME_${opt}} is not one of: ${allowed}")
    endif()
endforeach()

if(FLUME_OMP AND NOT FLUME_GPU_BACKEND STREQUAL "NONE")
    message(FATAL_ERROR "FLUME_OMP=ON cannot be combined with FLUME_GPU_BACKEND=${FLUME_GPU_BACKEND}")
endif()

if(FLUME_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
endif()

# --- AMReX (from the ERF submodule) ----------------------------------------
set(FLUME_AMREX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF/Submodules/AMReX)
if(NOT EXISTS ${FLUME_AMREX_DIR}/CMakeLists.txt)
    message(FATAL_ERROR
        "AMReX was not found in Submodules/ERF/Submodules/AMReX. Run:\n"
        "  git submodule update --init Submodules/ERF\n"
        "  git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX")
endif()

# The flood grid is two-dimensional. The build is three-dimensional with boxes
# one cell thick, so the same code can later link with ERF.
set(AMReX_SPACEDIM       3)
set(AMReX_PRECISION      ${FLUME_PRECISION})
set(AMReX_MPI            ${FLUME_MPI})
set(AMReX_OMP            ${FLUME_OMP})
set(AMReX_GPU_BACKEND    ${FLUME_GPU_BACKEND})
set(AMReX_FORTRAN        OFF)
set(AMReX_INSTALL        OFF)
set(AMReX_LINEAR_SOLVERS OFF)
set(AMReX_PARTICLES      OFF)
set(AMReX_AMRLEVEL       OFF)
set(AMReX_EB             OFF)
set(AMReX_FFT            OFF)
if(FLUME_TESTING)
    set(AMReX_TESTING    ON)
    set(AMReX_ASSERTIONS ON)
endif()

# SYSTEM: AMReX's headers are not held to FLUME's warning flags.
add_subdirectory(${FLUME_AMREX_DIR} ${CMAKE_BINARY_DIR}/amrex SYSTEM)
# For AMReXTargetHelpers (setup_target_for_cuda_compilation).
list(APPEND CMAKE_MODULE_PATH ${FLUME_AMREX_DIR}/Tools/CMake)

# MSVC: M_PI and friends need _USE_MATH_DEFINES; NOMINMAX keeps <windows.h>
# from defining min and max macros.
if(WIN32)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
endif()

add_subdirectory(Source)

# --- Tests ------------------------------------------------------------------
if(FLUME_ENABLE_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_subdirectory(Tests)
endif()
