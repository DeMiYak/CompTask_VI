cmake_minimum_required(VERSION 3.25)
project(CompTask_VI)

set(CMAKE_CXX_STANDARD 17)

include_directories(C:/MyFiles/ExtLib\ C++/boost_1_82_0)

add_executable(CompTask_VI main.cpp RPN.cpp SolveODE.cpp TaylorTask.cpp AdamsTask.cpp Runge-KuttaTask.cpp EulerTask.cpp)
