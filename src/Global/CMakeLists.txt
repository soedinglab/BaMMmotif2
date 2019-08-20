file(GLOB SOURCES *.cpp *.h ../getopt_pp/*)

add_executable (BaMMmotif ${SOURCES})

target_link_libraries (BaMMmotif LINK_PUBLIC init)
target_link_libraries (BaMMmotif LINK_PUBLIC refiner)
target_link_libraries (BaMMmotif LINK_PUBLIC evaluater)
target_link_libraries (BaMMmotif LINK_PUBLIC scanner)
target_link_libraries (BaMMmotif LINK_PUBLIC simulator)
target_link_libraries (BaMMmotif LINK_PUBLIC predictor)

add_library (global ${SOURCES})
