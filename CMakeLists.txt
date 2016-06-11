project(qahwa)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(includes)
add_subdirectory(sample)
add_subdirectory(src)
