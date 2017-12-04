file(GLOB SOURCES *.cpp *.h ../seq_generator/SeqGenerator.*)

add_executable (BaMMScan ${SOURCES})

target_link_libraries (BaMMScan LINK_PUBLIC init)

add_library (seq_scoring ${SOURCES})
