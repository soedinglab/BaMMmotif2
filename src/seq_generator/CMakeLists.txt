file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMSimu ${SOURCES})

target_link_libraries (BaMMSimu LINK_PUBLIC global)
target_link_libraries (BaMMSimu LINK_PUBLIC init)
target_link_libraries (BaMMSimu LINK_PUBLIC refiner)
target_link_libraries (BaMMSimu LINK_PUBLIC getopt_pp)

add_library (simulator ${SOURCES})
