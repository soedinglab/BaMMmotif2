file(GLOB SOURCES *.cpp *.h)
add_library (seq_generator ${SOURCES})

#target_include_directories (seq_generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})