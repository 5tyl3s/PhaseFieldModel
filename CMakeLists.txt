cmake_minimum_required(VERSION 3.10)
project(PhaseFieldModel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(C:/msys64/mingw64/include)

# Manual library paths for MSYS2 64-bit
link_directories(C:/msys64/mingw64/lib)

# Source files
set(SOURCES
    src/glVisualization.cpp
    src/gridField.cpp
    src/importConfig.cpp
    src/main.cpp
    src/model.cpp
    src/modelStartup.cpp
)

add_executable(PhaseFieldModel ${SOURCES})

# Link libraries directly
target_link_libraries(PhaseFieldModel PRIVATE
    glfw3
    glew32
    opengl32
)

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(PhaseFieldModel PRIVATE OpenMP::OpenMP_CXX)

