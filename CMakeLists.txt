cmake_minimum_required(VERSION 3.16)
project(matrix_utils CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(matrix src/matrix.cpp)
target_include_directories(matrix PUBLIC include)

add_executable(matrix_demo src/main.cpp)
target_link_libraries(matrix_demo matrix)

add_executable(matrix_tests tests/matrix_tests.cpp)
target_link_libraries(matrix_tests matrix)
