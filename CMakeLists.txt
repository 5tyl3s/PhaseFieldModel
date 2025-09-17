cmake_minimum_required(VERSION 3.10)
project(PhaseFieldModel)

# Set build type to Debug for better debugging information
set(CMAKE_BUILD_TYPE Debug)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# For MinGW OpenMP support
if(MINGW)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    # Use MinGW-specific console subsystem flag
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--subsystem,console")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Collect all source files
file(GLOB SOURCES 
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Print source files for debugging
message(STATUS "Source files: ${SOURCES}")

# Add executable
add_executable(PhaseFieldModel ${SOURCES})

# Add debug flags
target_compile_options(PhaseFieldModel PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3>
)

# Print executable location
get_target_property(EXEC_PATH PhaseFieldModel RUNTIME_OUTPUT_DIRECTORY)
message(STATUS "Executable will be created at: ${EXEC_PATH}")