cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -w") #Set warnings here, TODO
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")



set(CMAKE_PREFIX_PATH "/home/rasmus/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11")

find_package(Armadillo REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Python3 COMPONENTS Development REQUIRED)
find_package(mlpack REQUIRED)
find_package(Threads REQUIRED) #Parlay stuff
find_package(Parlay CONFIG REQUIRED) #Parlay stuff


# Ensure all code is compiled with -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})







include_directories(src)
add_executable(MyProject src/example.cpp
                src/kcentroids.cpp
                src/dc_hdbscan.cpp
                src/dc_dist.cpp
                src/dc_dist_bindings.cpp
                src/graph.cpp
                src/graph_hdb_c.cpp
                )


target_link_libraries(MyProject mlpack armadillo pybind11::module ${Python3_LIBRARIES} Parlay::parlay) #Parlay::parlay


# Add the static library - Needed to link for the pybind (here all cpp files needed that aren't the module itself should be added).
add_library(dcdist STATIC src/dc_dist.cpp src/dc_hdbscan.cpp src/graph.cpp src/graph_hdb_c.cpp
    )

# Add the Python module
pybind11_add_module(dctree src/dc_dist_bindings.cpp)

target_link_libraries(dctree PRIVATE dcdist ${MLPACK_LIBRARIES} ${Armadillo_LIBRARIES})
# Link the Python module to the static library

target_include_directories(${PROJECT_NAME} PRIVATE parallel_hdbscan/include)
target_include_directories(${PROJECT_NAME} PRIVATE parallel_hdbscan/src)



# Installation rules
install(TARGETS dctree LIBRARY DESTINATION "/usr/lib/python3/dist-packages")

