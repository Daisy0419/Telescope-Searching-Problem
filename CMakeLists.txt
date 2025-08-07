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

# # ---------- Gurobi ----------
# set(GUROBI_HOME "$ENV{GUROBI_HOME}")
# if(NOT GUROBI_HOME)
#     message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
# endif()
# set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
# set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
# set(GUROBI_LIB        "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
# set(GUROBI_CXX_LIB    "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a") # adjust if needed
# message(STATUS "Gurobi Library: ${GUROBI_LIB}")
# message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")
# set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})

# --- Gurobi paths (you likely already have these) ---
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
set(GUROBI_LIB_DIR "${GUROBI_HOME}/lib")
set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")

# Core solver .so
find_library(GUROBI_LIB NAMES gurobi120 HINTS "${GUROBI_LIB_DIR}" REQUIRED)

# Prefer shared C++ wrapper; fall back to old static archive
find_library(GUROBI_CXX_SHARED NAMES gurobi_c++ HINTS "${GUROBI_LIB_DIR}")
if (GUROBI_CXX_SHARED)
  set(GUROBI_CXX_LIB "${GUROBI_CXX_SHARED}")
  set(GRB_STATIC_WRAPPER OFF)
else()
  find_library(GUROBI_CXX_LIB NAMES gurobi_g++8.5 HINTS "${GUROBI_LIB_DIR}" REQUIRED)
  set(GRB_STATIC_WRAPPER ON)
endif()

# RPATH so runtime finds libgurobi120.so
list(APPEND CMAKE_BUILD_RPATH  "${GUROBI_LIB_DIR}")
list(APPEND CMAKE_INSTALL_RPATH "${GUROBI_LIB_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Common includes
set(PROJECT_INCLUDES
  ${CMAKE_SOURCE_DIR}/include
  ${GUROBI_INCLUDE_DIR}
  ${LEMON_INCLUDE_DIR}
)

# Link order MATTERS: C++ wrapper first, then core.
set(GRB_LINK_LIBS ${GUROBI_CXX_LIB} ${GUROBI_LIB})

# If using the static wrapper, defeat --as-needed (or group the libs)
if (GRB_STATIC_WRAPPER)
  add_link_options(-Wl,--no-as-needed)
  # Alternatively (either/or), you can group:
  # add_link_options(-Wl,--start-group) ; list(APPEND GRB_LINK_LIBS ${GUROBI_CXX_LIB} ${GUROBI_LIB}) ; add_link_options(-Wl,--end-group)
endif()

# # ---------- LEMON (adjust paths if needed) ----------
# set(LEMON_HOME "$ENV{LEMON_HOME}")
# if(NOT GUROBI_HOME)
#     message(LEMON_HOME "LEMON_HOME is not set. Please set it in your environment.")
# endif()
# set(LEMON_SOURCE_DIR "${GUROBI_HOME}")
# set(LEMON_BUILD_DIR  "${GUROBI_HOME}/build")

# ---------- LEMON ----------
find_path(LEMON_INCLUDE_DIR
  NAMES lemon/list_graph.h
  HINTS /usr/local/include /usr/include
  REQUIRED)

# Try common lib locations (lib and lib64)
find_library(LEMON_LIBRARY
  NAMES lemon
  HINTS /usr/local/lib /usr/local/lib64 /usr/lib /usr/lib/x86_64-linux-gnu
  )

# Fallback via pkg-config if the above didn’t find it
if (NOT LEMON_LIBRARY)
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    pkg_check_modules(LEMON QUIET lemon)
    if (LEMON_FOUND)
      set(LEMON_INCLUDE_DIR ${LEMON_INCLUDE_DIR} ${LEMON_INCLUDE_DIRS})
      set(LEMON_LIBRARY ${LEMON_LIBRARIES})
    endif()
  endif()
endif()

if (NOT LEMON_INCLUDE_DIR OR NOT LEMON_LIBRARY)
  message(FATAL_ERROR "LEMON not found. Checked headers at /usr/local/include and libs in /usr/local/lib{,64}.")
endif()

message(STATUS "LEMON include: ${LEMON_INCLUDE_DIR}")
message(STATUS "LEMON library: ${LEMON_LIBRARY}")


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

# ---------- Common libs ----------
set(PROJECT_LIBS
    OpenMP::OpenMP_CXX
    ${GUROBI_LIB}
    ${GUROBI_CXX_LIB}
    pthread m dl
    stdc++fs 
)

# ---------- Executable: ts (main.cpp) ----------
add_executable(ts ${MAIN} ${CORE_SOURCES})
# target_include_directories(ts PRIVATE ${PROJECT_INCLUDES})
# target_link_libraries(ts PRIVATE ${PROJECT_LIBS})

target_include_directories(ts PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(ts PRIVATE ${LEMON_LIBRARY} ${GRB_LINK_LIBS} OpenMP::OpenMP_CXX pthread m dl)


# ---------- Executable: op (main_custom.cpp) ----------
add_executable(op ${MAIN_CUSTOM} ${CORE_SOURCES})
# target_include_directories(op PRIVATE ${PROJECT_INCLUDES})
# target_link_libraries(op PRIVATE ${PROJECT_LIBS})

target_include_directories(op PRIVATE ${PROJECT_INCLUDES})
target_link_libraries(op PRIVATE ${LEMON_LIBRARY} ${GRB_LINK_LIBS} OpenMP::OpenMP_CXX pthread m dl)
