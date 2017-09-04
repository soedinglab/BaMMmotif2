file(GLOB SOURCES *.cpp *.h)
add_library (sequences ${SOURCES})

target_include_directories (sequences PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})