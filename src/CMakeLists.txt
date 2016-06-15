set(CPP_FILES
    Qahwa.cpp
    HistogramMaker.cpp
    ProbabilityDensityFunction.cpp
    UserDefinedFunctions.cpp
    UnbiasSolver.cpp
    )

set(QAHWA_INCLUDE_DIRECTORY "${PROJECT_SOURCE_DIR}/includes")
set(COFFEE_MILL_INCLUDE_DIRECTORY
    "${QAHWA_INCLUDE_DIRECTORY}/coffeemill/includes")
set(COFFEE_MILL_SOURCE_DIRECTORY
    "${QAHWA_INCLUDE_DIRECTORY}/coffeemill/src")

add_definitions("-O3 -Wall -std=c++11")
include_directories(${QAHWA_INCLUDE_DIRECTORY})
include_directories(${COFFEE_MILL_INCLUDE_DIRECTORY})
include_directories(${COFFEE_MILL_SOURCE_DIRECTORY})

add_executable(qahwa ${CPP_FILES})
target_link_libraries(qahwa
    coffeemill)
#     "${QAHWA_INCLUDE_DIRECTORY}/coffeemill/lib/libcoffeemill.a")
