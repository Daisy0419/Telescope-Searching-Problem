cmake_minimum_required(VERSION 3.11)
project(ts)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# # Enable profiling with gprof
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg") #-fno-inline
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

find_package(OpenMP REQUIRED)

# Gurobi config
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
if(NOT GUROBI_HOME)
    message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
endif()
set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
set(GUROBI_CXX_LIB "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")
message(STATUS "Gurobi Library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")
set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})

# LEMON config
set(LEMON_SOURCE_DIR "/home/research/w.yanwang/lemon-1.3.1")
set(LEMON_BUILD_DIR "/home/research/w.yanwang/lemon-1.3.1/build")

# Add subdirectory
add_subdirectory(Prize-Collecting-TSP-with-Budget)

# Main sources
set(SOURCES
    main.cpp
    helplers.cpp
    AntColony.cpp
    Greedy.cpp
    ReadData.cpp
    ILP_gurobi.cpp
    SimulatedAnnealing.cpp
    mst.cpp
    mstHoogeveen.cpp
    mstHoogeveenParallel.cpp
    Genetic.cpp
    Genetic_st.cpp
)

# executable
add_executable(ts ${SOURCES})

# set include directories
target_include_directories(ts PRIVATE
    ${GUROBI_INCLUDE_DIR}
    ${LEMON_SOURCE_DIR}
    ${LEMON_BUILD_DIR}
)

# link libraries
target_link_libraries(ts
    PRIVATE
        OpenMP::OpenMP_CXX
        ${GUROBI_LIB}
        ${GUROBI_CXX_LIB}
        pctsp_lib
        stdc++fs
)
