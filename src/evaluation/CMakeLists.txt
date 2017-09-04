file(GLOB SOURCES *.cpp *.h ../refinement/EM.* ../refinement/GibbsSampling.*)
add_executable (FDR ${SOURCES})

target_link_libraries (FDR LINK_PUBLIC init)
target_link_libraries (FDR LINK_PUBLIC seq_scoring)
target_link_libraries (FDR LINK_PUBLIC seq_generator)

add_library (fdr ${SOURCES})
