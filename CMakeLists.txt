cmake_minimum_required(VERSION 3.9.0)
project(QedFvCoef VERSION 0.1.0)

include(CTest)
include(GNUInstallDirs)
enable_testing()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_FLAGS "-W -Wall")
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -mtune=native -fno-math-errno")

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

add_subdirectory(lib)
add_subdirectory(cli)
