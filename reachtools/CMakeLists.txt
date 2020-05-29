cmake_minimum_required(VERSION 3.7)
project(ReachTools)

set(CMAKE_CXX_STANDARD 11)

add_executable(ReachTools main.cpp reachtools.cpp reachtools.h readgenome.cpp readgenome.h cxstring.cpp cxstring.h)
