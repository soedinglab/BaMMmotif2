cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(BaMMmotif)

set(BAMMMOTIF_MAJOR_VERSION 0)
set(BAMMMOTIF_MINOR_VERSION 1)
set(BAMMMOTIF_PATCH_VERSION 1)

set(BAMMMOTIF_VERSION ${BAMMMOTIF_MAJOR_VERSION}.${BAMMMOTIF_MINOR_VERSION}.${BAMMMOTIF_PATCH_VERSION})

SET (CMAKE_C_COMPILER /usr/bin/cc)
SET (CMKAE_CXX_COMPILER /usr/bin/c++)

set(CMAKE_CXX_FLAGS "-std=c++11 -DLOG_MAX_LEVEL=0 -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -O0 -g3 -pedantic -pedantic-errors -Wall -fmessage-length=0 -fno-strict-aliasing -Wconversion")

set(CMAKE_C_FLAGS "-D__GXX_EXPERIMENTAL_CXX0X__ -O0 -pedantic -pedantic-errors -Wall -fmessage-length=0 -Wconversion")
 
file(GLOB SOURCES "src/*.cpp" "src/*.h" "src/getopt_pp/*.h" "src/getopt_pp/*.cpp")
add_executable(BaMMmotif ${SOURCES})





