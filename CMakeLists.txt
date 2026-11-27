# ---------------------------------------------------------------------------
# FLUME -- top-level CMake build
#
#   git submodule update --init Submodules/ERF
#   git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX
#   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-rel -j 8
#   ctest --test-dir build-rel --output-on-failure
#
# AMReX is the one ERF pins: Submodules/ERF/Submodules/AMReX. FLUME has no
# AMReX pin of its own. The stand-alone build (the default) compiles no ERF
# code. With -DFLUME_WITH_ERF=ON, ERF is added as a subproject: it builds
# erf_exec and the one AMReX that FLUME and ERF both link.
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

option(FLUME_MPI          "Build with MPI"                                       ON)
option(FLUME_TESTING      "AMReX_TESTING + assertions (SNaN fills, as GPU CI)"   OFF)
option(FLUME_ENABLE_TESTS "Register the tests"                                   ON)
option(FLUME_WERROR       "Warnings in FLUME's own sources are errors"           OFF)
option(FLUME_WITH_ERF     "Also build ERF (erf_exec) from FLUME_ERF_DIR, on one AMReX" OFF)
set(FLUME_ERF_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF CACHE PATH
    "ERF source tree: AMReX is taken from its Submodules/AMReX; with FLUME_WITH_ERF it is also built")

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

# OpenMP is not supported: the kernels and tests are written and checked for
# MPI ranks and GPUs only. FLUME_OMP=ON, on the command line or left in a build
# cache, would otherwise be ignored without a word; OFF is harmless.
if(FLUME_OMP)
    message(FATAL_ERROR "FLUME_OMP=ON is no longer supported: FLUME builds without OpenMP. "
                        "Remove it from the command line or the cache (cmake -U FLUME_OMP).")
endif()

if(FLUME_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
endif()

# --- AMReX (from the ERF submodule), and ERF with FLUME_WITH_ERF -------------
set(FLUME_AMREX_DIR ${FLUME_ERF_DIR}/Submodules/AMReX)
if(NOT EXISTS ${FLUME_AMREX_DIR}/CMakeLists.txt)
    message(FATAL_ERROR
        "AMReX was not found in ${FLUME_AMREX_DIR} (FLUME_ERF_DIR=${FLUME_ERF_DIR}). Run:\n"
        "  git submodule update --init Submodules/ERF\n"
        "  git -C Submodules/ERF submodule update --init --depth 1 Submodules/AMReX")
endif()

if(FLUME_TESTING)
    set(AMReX_TESTING    ON)
    set(AMReX_ASSERTIONS ON)
endif()

if(FLUME_WITH_ERF)
    # ERF adds AMReX with its own options (CMake/SetAmrexOptions.cmake) and
    # FLUME links that AMReX. FLUME's options are the only place to set the
    # build's MPI, precision and GPU backend: an ERF cache entry that says
    # otherwise stops the configure instead of being overridden.
    if(NOT EXISTS ${FLUME_ERF_DIR}/CMake/BuildERFExe.cmake)
        message(FATAL_ERROR "FLUME_ERF_DIR: ${FLUME_ERF_DIR} is not an ERF source tree "
                            "(no CMake/BuildERFExe.cmake)")
    endif()
    if(NOT EXISTS ${FLUME_ERF_DIR}/Tests/ERFTestHelpers.cmake
       OR NOT EXISTS ${FLUME_ERF_DIR}/Source/FloodCoupling/ERF_FloodEngineRegistry.H)
        message(FATAL_ERROR "FLUME_ERF_DIR: ${FLUME_ERF_DIR} has no flood coupling or no "
                            "Tests/ERFTestHelpers.cmake; FLUME_WITH_ERF needs ERF from "
                            "erf-model/ERF#5375 (Nov 2026) or later")
    endif()
    set(_erf_cuda OFF)
    set(_erf_hip  OFF)
    set(_erf_sycl OFF)
    if(FLUME_GPU_BACKEND STREQUAL "CUDA")
        set(_erf_cuda ON)
    elseif(FLUME_GPU_BACKEND STREQUAL "HIP")
        set(_erf_hip ON)
    elseif(FLUME_GPU_BACKEND STREQUAL "SYCL")
        set(_erf_sycl ON)
    endif()
    # <ERF variable> <value> <FLUME option that sets it, or - for a fixed value>
    set(_erf_from_flume
        ERF_ENABLE_MPI    ${FLUME_MPI}       FLUME_MPI
        ERF_PRECISION     ${FLUME_PRECISION} FLUME_PRECISION
        ERF_ENABLE_CUDA   ${_erf_cuda}       FLUME_GPU_BACKEND
        ERF_ENABLE_HIP    ${_erf_hip}        FLUME_GPU_BACKEND
        ERF_ENABLE_SYCL   ${_erf_sycl}       FLUME_GPU_BACKEND
        ERF_ENABLE_OPENMP OFF                -
        ERF_ENABLE_FLOOD  ON                 -)
    list(LENGTH _erf_from_flume _n)
    math(EXPR _last "${_n} - 1")
    foreach(_i RANGE 0 ${_last} 3)
        math(EXPR _iv "${_i} + 1")
        math(EXPR _io "${_i} + 2")
        list(GET _erf_from_flume ${_i}  _var)
        list(GET _erf_from_flume ${_iv} _val)
        list(GET _erf_from_flume ${_io} _opt)
        set(_type BOOL)
        if(_var STREQUAL "ERF_PRECISION")
            set(_type STRING)
        endif()
        # FLUME's value and the cache value, with a boolean read as ON or OFF
        # (1, TRUE, YES...), so -DFLUME_MPI=1 and ERF_ENABLE_MPI=ON agree.
        if(_type STREQUAL "BOOL")
            if(_val)
                set(_val ON)
            else()
                set(_val OFF)
            endif()
        endif()
        set(_cached "$CACHE{${_var}}")
        if(_type STREQUAL "BOOL" AND DEFINED CACHE{${_var}})
            if(_cached)
                set(_cached ON)
            else()
                set(_cached OFF)
            endif()
        endif()
        # A cache entry that differs from FLUME's value and from the value
        # FLUME wrote last (FLUME_SET_<variable>) was set by the user.
        if(DEFINED CACHE{${_var}} AND NOT "${_cached}" STREQUAL "${_val}"
           AND NOT "${_cached}" STREQUAL "$CACHE{FLUME_SET_${_var}}")
            if(_opt STREQUAL "-")
                message(FATAL_ERROR "${_var}=$CACHE{${_var}}: FLUME_WITH_ERF always builds ERF with "
                                    "${_var}=${_val}. Remove it (cmake -U ${_var}).")
            endif()
            message(FATAL_ERROR "${_var}=$CACHE{${_var}} disagrees with FLUME's ${_opt}, "
                                "which sets it to ${_val} with FLUME_WITH_ERF. Set ${_opt} "
                                "instead and remove ${_var} (cmake -U ${_var}).")
        endif()
        # Written to the cache before ERF is added, so ERF's option() and
        # set(... CACHE ...) find it and keep it.
        set(${_var} ${_val} CACHE ${_type} "Set by FLUME (FLUME_WITH_ERF)" FORCE)
        set(FLUME_SET_${_var} ${_val} CACHE INTERNAL "")
    endforeach()

    # FLUME is compiled into erf_exec as a flood engine (erf.flood.engine =
    # flume, Source/ERFEngine). Engine directories the user lists in the
    # cache entry ERF_FLOOD_ENGINE_DIRS are kept and FLUME's is added to
    # them, for ERF only: the cache entry is not changed. It is created
    # first, so that ERF's set(... CACHE ...) leaves this list in place.
    set(ERF_FLOOD_ENGINE_DIRS "" CACHE STRING
        "Flood engine directories (each with ERFFloodEngine.cmake) to compile in; needs ERF_ENABLE_FLOOD")
    set(FLUME_ERF_ENGINE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Source/ERFEngine)
    set(ERF_FLOOD_ENGINE_DIRS $CACHE{ERF_FLOOD_ENGINE_DIRS})
    if(NOT FLUME_ERF_ENGINE_DIR IN_LIST ERF_FLOOD_ENGINE_DIRS)
        list(APPEND ERF_FLOOD_ENGINE_DIRS ${FLUME_ERF_ENGINE_DIR})
    endif()

    # SYSTEM, as for AMReX below: ERF's and AMReX's headers are not held to
    # FLUME's warning flags.
    add_subdirectory(${FLUME_ERF_DIR} ${CMAKE_BINARY_DIR}/ERF SYSTEM)
    if(NOT TARGET AMReX::amrex)
        message(FATAL_ERROR "ERF in ${FLUME_ERF_DIR} did not define AMReX::amrex; "
                            "FLUME_WITH_ERF needs ERF's internal AMReX (ERF_USE_INTERNAL_AMREX=ON)")
    endif()
else()
    # The flood grid is two-dimensional. The build is three-dimensional with
    # boxes one cell thick, so the same code links with ERF.
    set(AMReX_SPACEDIM       3)
    set(AMReX_PRECISION      ${FLUME_PRECISION})
    set(AMReX_MPI            ${FLUME_MPI})
    set(AMReX_OMP            OFF)
    set(AMReX_GPU_BACKEND    ${FLUME_GPU_BACKEND})
    set(AMReX_FORTRAN        OFF)
    set(AMReX_INSTALL        OFF)
    set(AMReX_LINEAR_SOLVERS OFF)
    set(AMReX_PARTICLES      OFF)
    set(AMReX_AMRLEVEL       OFF)
    set(AMReX_EB             OFF)
    set(AMReX_FFT            OFF)

    # SYSTEM: AMReX's headers are not held to FLUME's warning flags.
    add_subdirectory(${FLUME_AMREX_DIR} ${CMAKE_BINARY_DIR}/amrex SYSTEM)
endif()
# For AMReXTargetHelpers (setup_target_for_cuda_compilation).
list(APPEND CMAKE_MODULE_PATH ${FLUME_AMREX_DIR}/Tools/CMake)

# FLUME's executables in the top of the build tree. Set after ERF is added, so
# erf_exec stays where ERF's test helpers look for it.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

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
