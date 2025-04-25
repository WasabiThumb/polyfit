cmake_minimum_required(VERSION 3.29)
project(polyfit C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-O3 -lm")
set(CMAKE_CXX_FLAGS "-O3 -lm")

file(GLOB SRC_FILES src/*.c)
add_library(${PROJECT_NAME} STATIC ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC include)
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION 0.1.0 SOVERSION 0)
