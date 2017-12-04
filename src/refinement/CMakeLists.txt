file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMmotif ${SOURCES})

target_link_libraries (BaMMmotif LINK_PUBLIC getopt_pp)
target_link_libraries (BaMMmotif LINK_PUBLIC init)
target_link_libraries (BaMMmotif LINK_PUBLIC evaluation)
target_link_libraries (BaMMmotif LINK_PUBLIC seq_scoring)
target_link_libraries (BaMMmotif LINK_PUBLIC seq_generator)

add_library (bamm ${SOURCES})
