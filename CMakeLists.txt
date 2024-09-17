cmake_minimum_required(VERSION 3.26)
project(hopscotch)

set(CMAKE_CXX_STANDARD 20
)

add_executable(hopscotch main.cpp
        HopscotchHashSet.h
        tests.h
)
