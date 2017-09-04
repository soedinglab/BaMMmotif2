file(GLOB SOURCES *.cpp *.h)
add_library (init ${SOURCES})

target_include_directories (init PUBLIC ../sequences)