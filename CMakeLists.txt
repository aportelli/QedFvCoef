cmake_minimum_required(VERSION 3.11.0)
project(QedFvCoef VERSION 0.1.0)

include(CTest)
include(GNUInstallDirs)
enable_testing()

find_package(Git REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -march=native -mtune=native -fno-math-errno")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()



add_subdirectory(lib)
add_subdirectory(cli)
