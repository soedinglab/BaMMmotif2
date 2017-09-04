cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)

project (BaMMmotif-v2)

# link to Boost library
set (Boost_USE_STATIC_LIBS ON CACHE BOOL "use static libraries from Boost")
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost REQUIRED)
include_directories (${Boost_INCLUDEDIR})
link_libraries (${Boost_LIBRARIES})

# set the version of BaMM!motif
set (BAMMMOTIF_MAJOR_VERSION 0)
set (BAMMMOTIF_MINOR_VERSION 1)
set (BAMMMOTIF_PATCH_VERSION 1)
set (BAMMMOTIF_VERSION ${BAMMMOTIF_MAJOR_VERSION}.${BAMMMOTIF_MINOR_VERSION}.${BAMMMOTIF_PATCH_VERSION})

set (CMAKE_CXX_FLAGS
        "-D__cplusplus=201103L -D__GXX_EXPERIMENTAL_CXX0X__ -DLOG_MAX_LEVEL=0 -fmessage-length=0 -fno-strict-aliasing -g3 -O3 -pedantic -pedantic-errors -std=c++11 -Wall -Wconversion")

add_subdirectory (src/getopt_pp getopt_pp)
add_subdirectory (src/init init)
add_subdirectory (src/seq_generator seq_generator)
add_subdirectory (src/seq_scoring seq_scoring)
add_subdirectory (src/refinement bamm)
add_subdirectory (src/evaluation fdr)
