file(GLOB SOURCES *.cpp *.h)

add_executable (FDR ${SOURCES})

target_link_libraries (FDR LINK_PUBLIC global)
target_link_libraries (FDR LINK_PUBLIC init)
target_link_libraries (FDR LINK_PUBLIC refiner)
target_link_libraries (FDR LINK_PUBLIC scanner)
target_link_libraries (FDR LINK_PUBLIC simulator)
target_link_libraries (FDR LINK_PUBLIC getopt_pp)

add_library (evaluater ${SOURCES})
