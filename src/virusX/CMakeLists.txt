file(GLOB SOURCES *.cpp *.h)
add_executable (VirusX ${SOURCES})

target_link_libraries (VirusX LINK_PUBLIC getopt_pp)
target_link_libraries (VirusX LINK_PUBLIC shared)

