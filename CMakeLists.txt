# CMake version
cmake_minimum_required(VERSION 3.11)

# Project name
project(ts)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Set the build type to Release by default if not specified
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
# endif()

# # Set compiler flags
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

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
    ILP_scip.cpp
    ILP_gurobi.cpp
)

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
    ILP_scip.cpp
    ILP_gurobi.cpp
)

# Header files
set(HEADERS
    helplers.h
    AntColony.h
    Genetic.h
    Greedy.h
    kMST.h
    Christofides.h
    ReadData.h
    BranchBound.h
    ILP_scip.h
)

# Find Gurobi installation
find_path(GUROBI_INCLUDE_DIR NAMES gurobi_c++.h PATHS /Library/gurobi*/mac64/include)
find_library(GUROBI_LIBRARY NAMES gurobi_c++ gurobi PATHS /Library/gurobi*/mac64/lib)
set(GUROBI_INCLUDE_DIR "/Library/gurobi1201/macos_universal2/include")
# set(GUROBI_LIBRARY "/Library/gurobi1201/macos_universal2/lib/libgurobi_c++.a")
set(GUROBI_LIBRARY "/Library/gurobi1201/macos_universal2/lib/libgurobi120.dylib")
set(GUROBI_CPP_LIBRARY "/Library/gurobi1201/macos_universal2/lib/libgurobi_c++.a")


# set(GUROBI_LIBRARY "/Library/gurobi1201/macos_universal2/lib/libgurobi120.dylib")

message(STATUS "GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")
message(STATUS "GUROBI_LIBRARY: ${GUROBI_LIBRARY}")

if (NOT GUROBI_INCLUDE_DIR OR NOT GUROBI_LIBRARY)
    message(FATAL_ERROR "Gurobi not found. Make sure it is properly installed.")
endif()


# Find SCIP installation 
find_path(SCIP_INCLUDE_DIR NAMES scip/scip.h PATHS /usr/local/include)
find_library(SCIP_LIBRARY NAMES scip PATHS /usr/local/lib)
set(SCIP_INCLUDE_DIR "/usr/local/include")
set(SCIP_LIBRARY "/usr/local/lib/libscip.dylib")

message(STATUS "SCIP_INCLUDE_DIR: ${SCIP_INCLUDE_DIR}")
message(STATUS "SCIP_LIBRARY: ${SCIP_LIBRARY}")

if (NOT SCIP_INCLUDE_DIR OR NOT SCIP_LIBRARY)
    message(FATAL_ERROR "SCIP not found. Make sure it is properly installed.")
endif()


# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Specify include directories
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR} ${GUROBI_INCLUDE_DIR} ${SCIP_INCLUDE_DIR})

# Link Gurobi and SCIP libraries
# target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_LIBRARY} ${SCIP_LIBRARY})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GUROBI_LIBRARY} ${GUROBI_CPP_LIBRARY} ${SCIP_LIBRARY})
# target_link_libraries(${PROJECT_NAME} PRIVATE ${SCIP_LIBRARY})

# # CMake version
# cmake_minimum_required(VERSION 3.11)

# # Project name
# project(ts)

# # Set the C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Set the build type to Release by default if not specified
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
# endif()

# # Set compiler flags
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# # if(NOT CMAKE_BUILD_TYPE)
# #     message(STATUS "Setting build type to Debug as none was specified.")
# #     set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the build type: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
# # endif()

# # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# # Source files
# set(SOURCES
#     main.cpp
#     helplers.cpp
#     AntColony.cpp
#     Genetic.cpp
#     Greedy.cpp
#     kMST.cpp
#     Christofides.cpp
#     ReadData.cpp
#     BranchBound.cpp
#     ILP_scip.cpp
# )

# # Header files (not needed in add_executable)
# set(HEADERS
#     helplers.h
#     AntColony.h
#     Genetic.h
#     Greedy.h
#     kMST.h
#     Christofides.h
#     ReadData.h
#     BranchBound.h
#     ILP_scip.h
# )

# # Create executable
# add_executable(${PROJECT_NAME} ${SOURCES})

# # Specify include directories
# target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

