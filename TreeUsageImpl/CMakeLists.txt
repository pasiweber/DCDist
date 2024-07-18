cmake_minimum_required(VERSION 3.10)
project(MyProject)



find_package(Armadillo REQUIRED)


set(MLPACK_DIR external_libs/mlpack-4.4.0)
include_directories(${MLPACK_DIR}/include)
link_directories(${MLPACK_DIR}/lib)


include_directories(src)
add_executable(MyProject src/example.cpp
                src/kcentroids.cpp
                src/dc_hdbscan.cpp
                src/dc_dist.cpp
                )


target_link_libraries(MyProject mlpack armadillo)