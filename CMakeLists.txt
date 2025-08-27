cmake_minimum_required(VERSION 3.10)
project(FoilSimulation LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SFML 2 REQUIRED COMPONENTS graphics window system)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-window sfml-system)