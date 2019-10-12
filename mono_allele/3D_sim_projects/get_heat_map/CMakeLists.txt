cmake_minimum_required(VERSION 2.6)
project(get_heat_map)

add_executable(get_heat_map main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(get_heat_map PRIVATE -O0 -W -Wall -std=c++17)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(get_heat_map PRIVATE -O3 -W -Wall -std=c++17 -static)
target_link_libraries(get_heat_map PRIVATE stdc++fs -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

install(TARGETS get_heat_map RUNTIME DESTINATION bin)
