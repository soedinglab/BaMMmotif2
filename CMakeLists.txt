cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)

project (BaMMmotif-dev)

set (BAMMMOTIF_MAJOR_VERSION 0)
set (BAMMMOTIF_MINOR_VERSION 1)
set (BAMMMOTIF_PATCH_VERSION 1)

set (BAMMMOTIF_VERSION ${BAMMMOTIF_MAJOR_VERSION}.${BAMMMOTIF_MINOR_VERSION}.${BAMMMOTIF_PATCH_VERSION})

set (CMAKE_C_FLAGS "-D__GXX_EXPERIMENTAL_CXX0X__ -O3 -pedantic -pedantic-errors -Wall -fmessage-length=0 -Wconversion")

set (CMAKE_CXX_FLAGS "-D__cplusplus=201103L -D__GXX_EXPERIMENTAL_CXX0X__ -DLOG_MAX_LEVEL=0 -fmessage-length=0 -fno-strict-aliasing -g3 -O3 -pedantic -pedantic-errors -std=c++11 -Wall -Wconversion")

add_subdirectory (src/getopt_pp getopt_pp)
add_subdirectory (src/shared shared)
add_subdirectory (src/bamm BaMMmotif)
add_subdirectory (src/virusX VirusX)

