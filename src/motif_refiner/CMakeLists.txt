file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMRefiner ${SOURCES})

target_link_libraries (BaMMRefiner LINK_PUBLIC global)
target_link_libraries (BaMMRefiner LINK_PUBLIC init)
target_link_libraries (BaMMRefiner LINK_PUBLIC simulator)

add_library (refiner ${SOURCES})