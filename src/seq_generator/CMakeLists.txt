file(GLOB SOURCES *.cpp *.h ../refinement/EM.* )

add_executable (BaMMSimu ${SOURCES})

target_link_libraries (BaMMSimu LINK_PUBLIC init)

add_library (seq_generator ${SOURCES})
