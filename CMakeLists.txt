cmake_minimum_required(VERSION 3.29)
project(polyfit C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-O3 -lm -fPIC")
set(CMAKE_CXX_FLAGS "-O3 -lm -fPIC")

file(GLOB SRC_FILES src/*.c)
add_library(${PROJECT_NAME} SHARED ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC include)
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION 0.1.1 SOVERSION 0)
