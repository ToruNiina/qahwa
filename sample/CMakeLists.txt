set(QAHWA_INCLUDE_DIRECTORY "${PROJECT_SOURCE_DIR}/includes")
set(COFFEE_MILL_INCLUDE_DIRECTORY
    "${QAHWA_INCLUDE_DIRECTORY}/coffeemill/includes")
set(COFFEE_MILL_LIBRARY_INCLUDE_DIRECTORY
    "${QAHWA_INCLUDE_DIRECTORY}/coffeemill/src")

add_definitions("-O3 -Wall -std=c++11")
include_directories(${QAHWA_INCLUDE_DIRECTORY})
include_directories(${COFFEE_MILL_INCLUDE_DIRECTORY})
include_directories(${COFFEE_MILL_LIBRARY_INCLUDE_DIRECTORY})

add_executable(sample Sample.cpp)
set_target_properties(sample
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/sample"
    )
target_link_libraries(sample
    coffeemill)
