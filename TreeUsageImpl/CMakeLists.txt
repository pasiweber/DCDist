cmake_minimum_required(VERSION 3.10)
project(MyProject)
include_directories(src)
add_executable(MyProject src/example.cpp
                src/kcentroids.cpp
                src/dc_hdbscan.cpp
                src/dc_dist.cpp
                )