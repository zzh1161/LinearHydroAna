cmake_minimum_required(VERSION 3.10)

project(hex_cube)
set(CXX_STANDARD 17)

link_directories("/usr/lib/python3.8/config-3.8-x86_64-linux-gnu")
link_directories("/usr/loacl/lib")
include_directories("/usr/include/python3.8")

add_executable(main main.cpp)

target_link_libraries(main -lpython3.8 -lmshio)