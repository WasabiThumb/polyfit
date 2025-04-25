cmake_minimum_required(VERSION 3.29)
project(polyfit_test C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-O3 -lm")
set(CMAKE_CXX_FLAGS "-O3 -lm")

file(GLOB LIB_SRC_FILES ../src/*.c)
add_library(${PROJECT_NAME} main.c ${LIB_SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ../include)
