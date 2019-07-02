# CMakeLists.txt for event package. It creates a library with dictionary and a main program

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(production)

set (CMAKE_CXX_FLAGS "-std=c++14 -DVERBOSE=1")
set(BUILD_DIR ${PROJECT_SOURCE_DIR}/Build) # the Build area

# Set the output folder where your program will be created
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

####### ROOT ##########
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

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

## make install
install(FILES  ${CMAKE_CURRENT_BINARY_DIR}/Source/libproduction_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_PATH})

## make distclean ==> clean all build files and start from scratch
add_custom_target(distclean
    COMMAND find ${BUILD_DIR}/. ! -name 'build.sh' -type f -exec rm -f {} +
    COMMAND rm -f ${LIBRARY_OUTPUT_PATH}/lib*
    COMMAND rm -rf ${BUILD_DIR}/Source ${BUILD_DIR}/CMakeFiles
)


