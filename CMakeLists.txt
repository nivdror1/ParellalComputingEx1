cmake_minimum_required(VERSION 3.6)
project(ex1)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11++")

set(SOURCE_FILES main.cpp)
add_executable(ex1 ${SOURCE_FILES})