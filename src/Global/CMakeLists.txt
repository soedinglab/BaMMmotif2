file(GLOB SOURCES *.cpp *.h ../getopt_pp/*)

add_library (global ${SOURCES})
