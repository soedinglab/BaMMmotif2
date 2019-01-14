file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMmotif ${SOURCES})

target_link_libraries (BaMMmotif LINK_PUBLIC global)
target_link_libraries (BaMMmotif LINK_PUBLIC init)
target_link_libraries (BaMMmotif LINK_PUBLIC evaluater)
target_link_libraries (BaMMmotif LINK_PUBLIC scanner)
target_link_libraries (BaMMmotif LINK_PUBLIC simulator)

add_library (refiner ${SOURCES})