cmake_minimum_required(VERSION 2.8)
project(bin2mesh)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(ultimaille)
include_directories(ultimaille ultimaille/ext)


add_executable(bin2mesh main.cpp bin2mesh.cpp bin2mesh.h)
target_link_libraries(bin2mesh ultimaille)


