file(GLOB SOURCES *.cpp *.h)
add_library (seq_scoring ${SOURCES})

target_include_directories  (seq_scoring PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

