cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)

project (BaMMmotif-v2)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FindASan)

# link to Boost library
set (Boost_USE_STATIC_LIBS ON CACHE BOOL "use static libraries from Boost")
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost REQUIRED)
include_directories (${Boost_INCLUDE_DIR})
link_libraries (${Boost_LIBRARIES})

# link to Eigen library
set(EIGEN3_INCLUDE_DIR "~/opt/include/eigen3/")
if(NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( INSTALLATION_ERROR "Please install eigen3 using ./eigen3.downloader .")
endif()
include_directories (${EIGEN3_INCLUDE_DIR})

# set the version of BaMM!motif
set (BAMMMOTIF_MAJOR_VERSION 2)
set (BAMMMOTIF_MINOR_VERSION 0)
set (BAMMMOTIF_PATCH_VERSION 0)
set (BAMMMOTIF_VERSION ${BAMMMOTIF_MAJOR_VERSION}.${BAMMMOTIF_MINOR_VERSION}.${BAMMMOTIF_PATCH_VERSION})

set (CMAKE_CXX_FLAGS "-std=c++11 -Wall -g")

find_package(OpenMP)
if(OPENMP_FOUND OR OpenMP_FOUND OR OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP)
endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

install (DIRECTORY R/ ${CMAKE_BINARY_DIR}/bin/
        DESTINATION bin
        PATTERN *
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
        )

add_subdirectory (src/getopt_pp getopt_pp)
add_subdirectory (src/init init)
add_subdirectory (src/seq_generator simulator)
add_subdirectory (src/seq_scanner scanner)
add_subdirectory (src/refiner bamm)
add_subdirectory (src/evaluater evaluater)
add_subdirectory (src/Global global)
