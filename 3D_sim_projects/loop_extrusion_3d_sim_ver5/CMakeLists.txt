cmake_minimum_required(VERSION 2.6)
project(loop_extrusion_3d_sim_ver5)

add_executable(loop_extrusion_3d_sim_ver5 main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(loop_extrusion_3d_sim_ver5 PRIVATE -O0 -W -Wall -std=c++17)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(loop_extrusion_3d_sim_ver5 PRIVATE -O3 -W -Wall -std=c++17 -static)
target_link_libraries(loop_extrusion_3d_sim_ver5 PRIVATE stdc++fs -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

install(TARGETS loop_extrusion_3d_sim_ver5 RUNTIME DESTINATION bin)
