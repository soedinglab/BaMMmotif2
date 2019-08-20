file(GLOB SOURCES *.cpp *.h)

add_executable (SeqScanner ${SOURCES})

target_link_libraries (SeqScanner LINK_PUBLIC init)
target_link_libraries (SeqScanner LINK_PUBLIC global)
target_link_libraries (SeqScanner LINK_PUBLIC simulator)

add_library (scanner ${SOURCES})
