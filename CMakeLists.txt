cmake_minimum_required(VERSION 3.2)
project(mining)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")
set(SOURCE_FILES mining.cpp)
add_executable(mining ${SOURCE_FILES})