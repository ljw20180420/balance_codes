cmake_minimum_required(VERSION 2.6)
project(set_range_and_ctcf_ver3)

add_executable(set_range_and_ctcf_ver3 main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(set_range_and_ctcf_ver3 PRIVATE -O0 -W -Wall -std=c++17)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(set_range_and_ctcf_ver3 PRIVATE -O3 -W -Wall -std=c++17 -static)
target_link_libraries(set_range_and_ctcf_ver3 PRIVATE stdc++fs -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

install(TARGETS set_range_and_ctcf_ver3 RUNTIME DESTINATION bin)
