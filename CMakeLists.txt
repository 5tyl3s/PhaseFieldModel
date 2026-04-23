cmake_minimum_required(VERSION 3.10)
project(PhaseFieldModel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile-time option to enable/disable visualization
option(ENABLE_VISUALIZATION "Enable OpenGL visualization (requires GLFW, GLEW, OpenGL)" ON)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files (always included)
set(SOURCES
    src/gridField.cpp
    src/importConfig.cpp
    src/main.cpp
    src/model.cpp
    src/modelStartup.cpp
    src/helperFunctions.cpp
)

# Add visualization source only if enabled
if(ENABLE_VISUALIZATION)
    list(APPEND SOURCES src/glVisualization.cpp)
endif()

add_executable(PhaseFieldModel ${SOURCES})

# Add compile definition for visualization toggle
if(ENABLE_VISUALIZATION)
    target_compile_definitions(PhaseFieldModel PRIVATE ENABLE_VISUALIZATION=1)
else()
    target_compile_definitions(PhaseFieldModel PRIVATE ENABLE_VISUALIZATION=0)
endif()

# OpenMP (required)
find_package(OpenMP REQUIRED)
target_link_libraries(PhaseFieldModel PRIVATE OpenMP::OpenMP_CXX)

# Optional: OpenGL/GLFW for visualization (only if enabled)
# Try multiple approaches to find these libraries

if(ENABLE_VISUALIZATION)
    # First, try pkg-config (works on many systems)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(GLFW glfw3 QUIET)
        pkg_check_modules(GLEW glew QUIET)
    endif()

    # If pkg-config didn't find them, try CMake find_package
    if(NOT GLFW_FOUND)
        find_package(GLFW3 QUIET)
        if(GLFW3_FOUND)
            set(GLFW_LIBRARIES ${GLFW3_LIBRARIES})
            set(GLFW_INCLUDE_DIRS ${GLFW3_INCLUDE_DIR})
        endif()
    endif()

    if(NOT GLEW_FOUND)
        find_package(GLEW QUIET)
        if(GLEW_FOUND)
            set(GLEW_LIBRARIES ${GLEW_LIBRARIES})
            set(GLEW_INCLUDE_DIRS ${GLEW_INCLUDE_DIRS})
        endif()
    endif()

    find_package(OpenGL QUIET)

    # Check if we found the visualization libraries
    if(GLFW_FOUND AND GLEW_FOUND AND OPENGL_FOUND)
        target_link_libraries(PhaseFieldModel PRIVATE ${GLFW_LIBRARIES} ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})
        if(GLFW_INCLUDE_DIRS)
            target_include_directories(PhaseFieldModel PRIVATE ${GLFW_INCLUDE_DIRS})
        endif()
        if(GLEW_INCLUDE_DIRS)
            target_include_directories(PhaseFieldModel PRIVATE ${GLEW_INCLUDE_DIRS})
        endif()
        message(STATUS "Visualization libraries found via pkg-config/CMake - OpenGL support enabled")
    else()
        # Fallback: try manual MSYS2 paths
        message(STATUS "Trying MSYS2 fallback paths...")
        set(MSYS2_LIB_PATHS 
            "C:/msys64/mingw64/lib"
            "C:/msys64/mingw32/lib"
            "C:/mingw64/lib"
            "C:/mingw32/lib"
        )
        
        set(MSYS2_INCLUDE_PATHS
            "C:/msys64/mingw64/include"
            "C:/msys64/mingw32/include"
            "C:/mingw64/include"
            "C:/mingw32/include"
        )
        
        find_library(GLFW3_LIB NAMES glfw3 glfw PATHS ${MSYS2_LIB_PATHS} NO_DEFAULT_PATH)
        find_library(GLEW_LIB NAMES glew32 glew PATHS ${MSYS2_LIB_PATHS} NO_DEFAULT_PATH)
        find_library(OPENGL_LIB NAMES opengl32 opengl PATHS ${MSYS2_LIB_PATHS} NO_DEFAULT_PATH)
        
        if(GLFW3_LIB AND GLEW_LIB AND OPENGL_LIB)
            target_link_libraries(PhaseFieldModel PRIVATE ${GLFW3_LIB} ${GLEW_LIB} ${OPENGL_LIB})
            
            # Add include directories from MSYS2
            foreach(INCLUDE_PATH ${MSYS2_INCLUDE_PATHS})
                if(EXISTS ${INCLUDE_PATH})
                    target_include_directories(PhaseFieldModel PRIVATE ${INCLUDE_PATH})
                endif()
            endforeach()
            
            message(STATUS "Visualization libraries found via MSYS2 fallback - OpenGL support enabled")
        else()
            message(WARNING "Visualization libraries NOT found")
            message(WARNING "  - GLFW3_LIB: ${GLFW3_LIB}")
            message(WARNING "  - GLEW_LIB: ${GLEW_LIB}")
            message(WARNING "  - OPENGL_LIB: ${OPENGL_LIB}")
            message(WARNING "Build will succeed but set enableVisualization=false in config.json")
        endif()
    endif()
else()
    message(STATUS "Visualization disabled at compile time")
endif()

