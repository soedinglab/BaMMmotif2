file(GLOB SOURCES *.cpp *.h ../refinement/utils.h)

add_executable (extractProbs ${SOURCES})

target_link_libraries (extractProbs LINK_PUBLIC init)

