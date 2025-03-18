cmake_minimum_required(VERSION 3.11)
project(ts)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Set Gurobi paths manually
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")

if(NOT GUROBI_HOME)
    message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
endif()

# Explicitly define the Gurobi libraries
set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
set(GUROBI_CXX_LIB "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")

# Debug messages
message(STATUS "Gurobi Library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")

# Source files
set(SOURCES
    main.cpp
    helplers.cpp
    AntColony.cpp
    Genetic.cpp
    Greedy.cpp
    kMST.cpp
    Christofides.cpp
    ReadData.cpp
    BranchBound.cpp
    ILP_gurobi.cpp
)

# Create the executable
add_executable(ts ${SOURCES})

# Include Gurobi headers
target_include_directories(ts PRIVATE ${GUROBI_INCLUDE_DIR})

# Link Gurobi libraries
# target_link_libraries(ts OpenMP::OpenMP_CXX ${GUROBI_LIB} ${GUROBI_CXX_LIB})
target_link_libraries(ts OpenMP::OpenMP_CXX ${GUROBI_LIB} ${GUROBI_CXX_LIB} stdc++fs)


# Ensure Gurobi runtime libraries are found
set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})
