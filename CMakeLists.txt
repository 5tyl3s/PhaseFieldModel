cmake_minimum_required(VERSION 3.10)
project(PhaseFieldModel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SOURCES "src/*.cpp")

add_executable(PhaseFieldModel 
    src/main.cpp 
    src/modelStartup.cpp
    ${SOURCES})