cmake_minimum_required(VERSION 3.11)
project(ts)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Profiling 
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# OpenMP
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

# CPLEX config
set(CPLEX_ROOT "/project/engineering/cplex22/cplex")
set(CPLEX_LIB_DIR "${CPLEX_ROOT}/lib/x86-64_linux/static_pic")
set(CONCERT_LIB_DIR "/project/engineering/cplex22/concert/lib/x86-64_linux/static_pic")
set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT}/include")
set(CPLEX_CONCERT_INCLUDE_DIR "/project/engineering/cplex22/concert/include")
set(ILOCLEX_LIB "${CPLEX_LIB_DIR}/libilocplex.a")
set(CPLEX_LIB "${CPLEX_LIB_DIR}/libcplex.a")
set(CONCERT_LIB "${CONCERT_LIB_DIR}/libconcert.a")

# LEMON config
set(LEMON_SOURCE_DIR "/home/research/w.yanwang/lemon-1.3.1")
set(LEMON_BUILD_DIR "/home/research/w.yanwang/lemon-1.3.1/build")

# Collect all source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Add executable
add_executable(ts ${SOURCES})

# Include directories
target_include_directories(ts PRIVATE
    include
    ${GUROBI_INCLUDE_DIR}
    ${CPLEX_INCLUDE_DIR}
    ${CPLEX_CONCERT_INCLUDE_DIR}
    ${LEMON_SOURCE_DIR}
    ${LEMON_BUILD_DIR}
)

# Link libraries
target_link_libraries(ts
    PRIVATE
        OpenMP::OpenMP_CXX
        ${GUROBI_LIB}
        ${GUROBI_CXX_LIB}
        ${ILOCLEX_LIB}
        ${CONCERT_LIB}
        ${CPLEX_LIB}
        pthread m dl
        stdc++fs
)
