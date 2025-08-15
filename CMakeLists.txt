cmake_minimum_required(VERSION 3.11)
project(op)

# ---------- C++17 ----------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---------- OpenMP ----------
find_package(OpenMP REQUIRED)

# ---------- Gurobi ----------
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
if(NOT GUROBI_HOME)
    message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
endif()
set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
find_library(GUROBI_LIB NAMES gurobi120 HINTS "${GUROBI_LIBRARY_DIR}" REQUIRED)
find_library(GUROBI_CXX_LIB NAMES gurobi_g++8.5 HINTS "${GUROBI_LIBRARY_DIR}" REQUIRED)
message(STATUS "Gurobi Library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")
set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})

# ---------- LEMON (adjust paths) ----------
#set(LEMON_SOURCE_DIR "/usr/local")
#set(LEMON_BUILD_DIR  "/usr/local")

# ---------- Sources ----------
# Grab all .cpp files
file(GLOB_RECURSE ALL_SOURCES "src/*.cpp")

# Explicitly identify the two mains
set(MAIN          "${CMAKE_SOURCE_DIR}/src/main.cpp")
set(MAIN_CUSTOM   "${CMAKE_SOURCE_DIR}/src/main_custom.cpp")

# Guard against missing files
if(NOT EXISTS "${MAIN}")
    message(FATAL_ERROR "Expected main at ${MAIN} not found.")
endif()
if(NOT EXISTS "${MAIN_CUSTOM}")
    message(FATAL_ERROR "Expected main at ${MAIN_CUSTOM} not found.")
endif()

# Core sources = all sources minus the two mains
set(CORE_SOURCES ${ALL_SOURCES})
list(REMOVE_ITEM CORE_SOURCES "${MAIN}" "${MAIN_CUSTOM}")

# ---------- Common include dirs ----------
set(PROJECT_INCLUDES
    ${CMAKE_SOURCE_DIR}/include
    ${GUROBI_INCLUDE_DIR}
    ${LEMON_SOURCE_DIR}
    ${LEMON_BUILD_DIR}
)

set(GRB_LINK_LIBS ${GUROBI_CXX_LIB} ${GUROBI_LIB})

# ---------- Common libs ----------
set(PROJECT_LIBS
    OpenMP::OpenMP_CXX
    ${GRB_LINK_LIBS}
    pthread m dl
    stdc++fs 
)

# ---------- Executable: ts (main.cpp) ----------
add_executable(ts ${MAIN} ${CORE_SOURCES})
target_include_directories(ts PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(ts PRIVATE ${PROJECT_LIBS})

# ---------- Executable: op (main_custom.cpp) ----------
add_executable(op ${MAIN_CUSTOM} ${CORE_SOURCES})
target_include_directories(op PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(op PRIVATE ${PROJECT_LIBS})