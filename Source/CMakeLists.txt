# The solver is a library so the executable and later unit tests link the
# same compiled code.
add_library(flume_lib STATIC
    FLUME_Advance.cpp
    FLUME_Boundary.cpp
    FLUME_Buildings.cpp
    FLUME_Flood.cpp
    FLUME_Friction.cpp
    FLUME_IO.cpp
    FLUME_Init.cpp
    FLUME_Inputs.cpp
    FLUME_Raster.cpp
    FLUME_Source.cpp
    FLUME_TimeStep.cpp)

target_include_directories(flume_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(flume_lib PUBLIC AMReX::amrex)

add_executable(flume main.cpp)
target_link_libraries(flume PRIVATE flume_lib)

# Warning flags for FLUME's own sources only (AMReX is a SYSTEM subdirectory).
foreach(tgt flume_lib flume)
    if(MSVC)
        target_compile_options(${tgt} PRIVATE /W4)
        if(FLUME_WERROR)
            target_compile_options(${tgt} PRIVATE /WX)
        endif()
    elseif(FLUME_GPU_BACKEND STREQUAL "NONE")
        target_compile_options(${tgt} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
        if(FLUME_WERROR)
            target_compile_options(${tgt} PRIVATE -Werror)
        endif()
    endif()
endforeach()

if(FLUME_GPU_BACKEND STREQUAL "CUDA")
    include(AMReXTargetHelpers)
    setup_target_for_cuda_compilation(flume_lib)
    setup_target_for_cuda_compilation(flume)
endif()
