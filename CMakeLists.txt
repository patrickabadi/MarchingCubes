cmake_minimum_required(VERSION 3.10)

if(POLICY CMP0074)
  # TODO: update *_ROOT variables to be PCL_*_ROOT or equivalent.
  # CMP0074 directly affects how Find* modules work and *_ROOT variables.  Since
  # this is a config file that will be consumed by parent projects with (likely)
  # NEW behavior, we need to push a policy stack.
  cmake_policy(PUSH)
  cmake_policy(SET CMP0074 OLD)
endif()

project(MarchingCubes VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FLANN_DIR $ENV{PCL_ROOT}share/flann)
set(lz4_DIR $ENV{PCL_ROOT}share/lz4)
set(PCL_DIR $ENV{PCL_ROOT}share/pcl)

find_package(OpenGL REQUIRED)

find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
  
add_subdirectory(MarchingCubesTest)

if(POLICY CMP0074)
  cmake_policy(POP)
endif()