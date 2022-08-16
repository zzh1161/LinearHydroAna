cmake_minimum_required(VERSION 3.10)

project(hex_smallCube)
set(CXX_STANDARD 17)

link_directories("/usr/loacl/lib")
include_directories("/usr/include/eigen3")
include_directories("/usr/local/include")

add_executable(main main.cpp)

target_link_libraries(main -lmshio)