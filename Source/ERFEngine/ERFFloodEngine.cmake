# FLUME as a flood engine of ERF (erf.flood.engine = flume). FLUME's
# CMakeLists.txt adds this directory to ERF_FLOOD_ENGINE_DIRS with
# FLUME_WITH_ERF; ERF includes this file and calls erf_flood_engine_add with
# its library. The engine links FLUME's solver library (flume_lib), which
# links the same AMReX as ERF.
function(erf_flood_engine_add target)
    target_sources(${target} PRIVATE ${ERF_FLOOD_ENGINE_DIR}/FLUME_ERFFloodEngine.cpp)
    target_link_libraries(${target} PUBLIC flume_lib)
endfunction()
