cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") #Set warnings here, TODO
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")



set(CMAKE_PREFIX_PATH "/home/rasmus/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11")

find_package(Armadillo REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Python3 COMPONENTS Development REQUIRED)

# Ensure all code is compiled with -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})


set(MLPACK_DIR external_libs/mlpack-4.4.0)
include_directories(${MLPACK_DIR}/include)
link_directories(${MLPACK_DIR}/lib)


find_package(Threads REQUIRED) #Parlay stuff
find_package(Parlay CONFIG REQUIRED) #Parlay stuff


include_directories(src)
add_executable(MyProject src/example.cpp
                src/kcentroids.cpp
                src/dc_hdbscan.cpp
                src/dc_dist.cpp
                src/dc_dist_bindings.cpp
                parallel_hdbscan/src/hdbscan_arma.cpp
                parallel_hdbscan/src/hdbscan.cpp
                parallel_hdbscan/src/dendrogram.cpp

                )


target_link_libraries(MyProject mlpack armadillo pybind11::module ${Python3_LIBRARIES} Parlay::parlay) #Parlay::parlay


# Add the static library (if needed) yeah
add_library(dcdist STATIC src/dc_dist.cpp)

# Add the Python module
pybind11_add_module(dctree src/dc_dist_bindings.cpp)

# Link the Python module to the static library
target_link_libraries(dctree PRIVATE dcdist ${MLPACK_LIBRARIES} ${Armadillo_LIBRARIES})

# Installation rules
install(TARGETS dctree LIBRARY DESTINATION "/usr/lib/python3/dist-packages")

