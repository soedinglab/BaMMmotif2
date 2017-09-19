file(GLOB SOURCES *.cpp *.h)
add_executable (BaMMScan ${SOURCES})

target_link_libraries (BaMMScan LINK_PUBLIC init)

add_library (seq_scoring ${SOURCES})
