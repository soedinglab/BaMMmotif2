file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMScan ${SOURCES})

target_link_libraries (BaMMScan LINK_PUBLIC init)
target_link_libraries (BaMMScan LINK_PUBLIC global)
target_link_libraries (BaMMScan LINK_PUBLIC simulator)

add_library (scanner ${SOURCES})
