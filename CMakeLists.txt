cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)

project (BaMMmotif-v2)

# link to Boost
set(Boost_USE_STATIC_LIBS ON CACHE BOOL "use static libraries from Boost") 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost REQUIRED) 
include_directories(${Boost_INCLUDEDIR})
link_libraries(${Boost_LIBRARIES})

set (BAMMMOTIF_MAJOR_VERSION 0)
set (BAMMMOTIF_MINOR_VERSION 1)
set (BAMMMOTIF_PATCH_VERSION 1)

set (BAMMMOTIF_VERSION ${BAMMMOTIF_MAJOR_VERSION}.${BAMMMOTIF_MINOR_VERSION}.${BAMMMOTIF_PATCH_VERSION})

set (CMAKE_C_FLAGS "-D__GXX_EXPERIMENTAL_CXX0X__ -O3 -pedantic -pedantic-errors -Wall -fmessage-length=0 -Wconversion")

set (CMAKE_CXX_FLAGS "-D__cplusplus=201103L -D__GXX_EXPERIMENTAL_CXX0X__ -DLOG_MAX_LEVEL=0 -fmessage-length=0 -fno-strict-aliasing -g3 -O3 -pedantic -pedantic-errors -std=c++11 -Wall -Wconversion -static-libstdc++")

file(GLOB SOURCES "src/*.h" "src/*.cpp" "src/getopt_pp/*.h" "src/getopt_pp/*.cpp")

add_executable (BaMMmotif ${SOURCES})

