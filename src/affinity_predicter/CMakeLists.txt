file(GLOB SOURCES *.cpp *.h)

add_executable (KmerPred ${SOURCES})

target_link_libraries (KmerPred LINK_PUBLIC init)
target_link_libraries (KmerPred LINK_PUBLIC global)

add_library (predictor ${SOURCES})
