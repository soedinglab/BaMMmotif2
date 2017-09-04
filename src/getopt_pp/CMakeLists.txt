file(GLOB SOURCES *.cpp *.h)
add_library (getopt_pp ${SOURCES})

target_include_directories  (getopt_pp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

