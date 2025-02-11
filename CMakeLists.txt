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
)

# Find SCIP library
find_library(SCIP_LIBRARY scip PATHS /usr/local/scip/lib REQUIRED)

# Include SCIP headers
include_directories(/usr/local/scip/include)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link only the SCIP library (without LPI)
target_link_libraries(${PROJECT_NAME} ${SCIP_LIBRARY})



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

