file(GLOB SOURCES *.cpp *.h)

add_executable (BaMMPred ${SOURCES})

target_link_libraries (BaMMPred LINK_PUBLIC init)
target_link_libraries (BaMMPred LINK_PUBLIC global)

add_library (predictor ${SOURCES})
