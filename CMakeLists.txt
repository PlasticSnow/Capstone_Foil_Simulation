cmake_minimum_required(VERSION 3.10)
project(FoilSimulation LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/sfml@2/lib/cmake/SFML/")

find_package(SFML 2 REQUIRED COMPONENTS graphics window system)

# Automatically gather all .cpp files in src/
file(GLOB SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-window sfml-system)

