# CMakeLists.txt for event package. It creates a library with dictionary and a main program

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(production)

set (CMAKE_CXX_FLAGS "-std=c++14 -DVERBOSE=1")

# Set the output folder where your program will be created
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

####### ROOT ##########
option(BUILD_WITH_ROOT
  "Enable ROOT"
  ON)
if(BUILD_WITH_ROOT)
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    find_package(ROOT REQUIRED COMPONENTS RIO Net)
    include(${ROOT_USE_FILE})
endif()

##########  CUDA ##########
option(BUILD_WITH_CUDA
  "Enable CUDA"
  OFF)

if(BUILD_WITH_CUDA)
  find_package(CUDA REQUIRED)
  include_directories(BEFORE SYSTEM "${CUDA_INCLUDE_DIRS}")
endif()

### compile options (note: not sure if this is necessary) ####
include_directories(BEFORE "${PROJECT_SOURCE_DIR}")


##### Subdirectories ######
##  need to be after above user options #######

add_subdirectory(Source)
