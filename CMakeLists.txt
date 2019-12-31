cmake_minimum_required(VERSION 3.15)
project(RTBP_PeriodicOrbitsL3 C)

# Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O")

# Source files
set(SOURCE_FILES main.c src/equPoints.c myHeader.h src/orbitIntegration.c src/IOmodule.c src/poincareSections.c src/symmetricPeriodicOrbits.c src/rtbpEqs.c src/rtbp.c src/rtbp_variationals.c src/taylor.h src/stability.c src/hillsRegions.c)
add_executable(RTBP_PeriodicOrbitsL3 ${SOURCE_FILES})

# External libraries:
# 1) math.h, linked using flag -lm
# 2) gsl (gnu scientific library), linked using -lgsl -lgslblas
target_link_libraries(RTBP_PeriodicOrbitsL3 m gsl gslcblas)