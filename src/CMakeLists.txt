file(GLOB SOURCES *.cpp *.h)
add_library (shared ${SOURCES})

target_include_directories (shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

