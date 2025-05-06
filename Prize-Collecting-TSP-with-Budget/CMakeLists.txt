# # Define source files
# set(PCTSP_SOURCES
#     Graph.cpp
#     PD.cpp
#     ReadFile.cpp
#     Subroutine.cpp
#     Subset.cpp
#     christofides.cpp
# )

# # Optional: define headers for clarity (not required by CMake)
# set(PCTSP_HEADERS
#     Graph.h
#     PD.h
#     ReadFile.h
#     Subroutine.h
#     Subset.h
#     christofides.h
# )

# # Export source variable to parent scope
# set(PCTSP_SOURCES ${PCTSP_SOURCES} PARENT_SCOPE)


# Source files in this subdirectory
set(PCTSP_SOURCES
    Graph.cpp
    PD.cpp
    ReadFile.cpp
    Subroutine.cpp
    Subset.cpp
    bpc_tsp.cpp
    christofides.cpp
)

add_library(pctsp_lib ${PCTSP_SOURCES})
target_include_directories(pctsp_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
