file(GLOB SOURCES *.cpp *.h)
add_executable (BaMMmotif ${SOURCES})

target_link_libraries (BaMMmotif LINK_PUBLIC getopt_pp)
target_link_libraries (BaMMmotif LINK_PUBLIC shared)

