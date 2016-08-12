file(GLOB SOURCES *Building.cpp *Building.h)
add_executable (build ${SOURCES})

target_link_libraries (build LINK_PUBLIC getopt_pp)
target_link_libraries (build LINK_PUBLIC shared)

file(GLOB SOURCES *Prediction.cpp *Prediction.h *Interactions.*)
add_executable (predict ${SOURCES})

target_link_libraries (predict LINK_PUBLIC getopt_pp)
target_link_libraries (predict LINK_PUBLIC shared)

