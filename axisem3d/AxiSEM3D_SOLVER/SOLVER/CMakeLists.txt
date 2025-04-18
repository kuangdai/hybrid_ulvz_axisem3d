#
#  CMakeLists.txt
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 15/8/19.
#  Copyright © 2019 Kuangdai Leng. All rights reserved.
#

#  CMakeLists.txt

########################################################################
####################### Build Settings (EDIT ME) #######################
########################################################################

# compiler suite
set(CMAKE_C_COMPILER           mpicc)
set(CMAKE_CXX_COMPILER         mpicxx)
set(CMAKE_Fortran_COMPILER     gfortran)

# preferred installation prefix of external packages
# ** the roots can also be set by environment variables, e.g.,
#    $ export EIGEN3_ROOT=$HOME/Developer/eigen3
# ** by default, the following CMake variables are initialized by
#    the environment variables; if specified otherwise, they will
#    have priority over the environment variables
set(EIGEN3_ROOT                $ENV{EIGEN3_ROOT})
set(BOOST_ROOT                 $ENV{BOOST_ROOT})
set(FFTW_ROOT                  $ENV{FFTW_ROOT})
set(METIS_ROOT                 $ENV{METIS_ROOT})
set(NETCDF_ROOT                $ENV{NETCDF_ROOT})
set(HDF5_ROOT                  $ENV{HDF5_ROOT})


# polynomial order of spectral elements (from 1 to 8)
set(NPOL 4)

# solver precision
# ** using double precision increases memory usage but barely
#    affects solver performance
set(USE_DOUBLE false)

# memory saving mode
# ** minimize memory usage with some performance hit
set(MEMORY_SAVING_MODE false)

# serial build
# ** build AxiSEM3D without MPI
set(SERIAL_BUILD false)

# use parallel NetCDF for output
# ** Par-NetCDF is supported but not required
set(USE_PARALLEL_NETCDF false)

# force link to HDF5 libraries
# ** set to "true" if linking fails with missing HDF5 symbols
set(FORCE_HDF5_LINK true)

# additional libraries to link to
# ** add libs if linking fails with missing symbols, e.g., "-lcurl"
set(ADDITIONAL_LIBS "-lcurl")

########################################################################
####################### Build Settings (EDIT ME) #######################
########################################################################



########### NOTE: Users are less likely to edit lines below ############



################# cmake setup #################
cmake_minimum_required(VERSION 3.0.0)
project(AxiSEM3D)
enable_language(Fortran)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# AxiSEM3D requires C++17
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# use "Release" as default build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# CMP0074 issue since CMake 3.12
if(POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()


################# macros passed to AxiSEM3D #################
# version
add_definitions(-D_VERSION="1.2")

# source dir
add_definitions(-D_PROJ_DIR=\"${PROJECT_SOURCE_DIR}\")

# nPol
add_definitions(-D_NPOL=${NPOL})

# double
if(USE_DOUBLE)
add_definitions(-D_USE_DOUBLE)
endif()

# memory
if(MEMORY_SAVING_MODE)
add_definitions(-D_SAVE_MEMORY)
endif()

# serial build
if(SERIAL_BUILD)
add_definitions(-D_SERIAL_BUILD)
set(USE_PARALLEL_NETCDF false)
endif()

# parallel NetCDF
if(USE_PARALLEL_NETCDF)
add_definitions(-D_USE_PARALLEL_NETCDF)
endif()


################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# mpi
if(NOT SERIAL_BUILD)
find_package(MPI COMPONENTS C REQUIRED)
include_directories(${MPI_C_INCLUDE_DIRS})
endif()

# eigen
# FindEIGEN3.cmake only uses ENV in HINTS
set(ENV{EIGEN3_ROOT} "${EIGEN3_ROOT}")
find_package(EIGEN3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# fftw
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

# metis
# FindMETIS.cmake uses METIS_DIR instead of METIS_ROOT in HINTS
set(METIS_DIR "${METIS_ROOT}")
set(ENV{METIS_DIR} "$ENV{METIS_ROOT}")
find_package(METIS REQUIRED QUIET)
include_directories(${METIS_INCLUDE_DIRS})

# NetCDF
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})

# hdf5
if(FORCE_HDF5_LINK)
find_package(HDF5 COMPONENTS C COMPONENTS HL REQUIRED)
endif()


################# include paths #################
# find all directories containing .hpp
# 1) scan *.hpp files in src/
file(GLOB_RECURSE HPP_PATH_FILES src/*.hpp)
# 2) separate path and filename
set(HPP_PATHS "")
foreach(HPP_PATH_FILE ${HPP_PATH_FILES})
get_filename_component(HPP_PATH ${HPP_PATH_FILE} PATH)
set(HPP_PATHS ${HPP_PATHS} ${HPP_PATH})
endforeach()
# 3) remove duplicates
list(REMOVE_DUPLICATES HPP_PATHS)

# include all directories containing .hpp
include_directories(
${HPP_PATHS}
# local externals
external
)


################# source files #################
# find *.cpp
file(GLOB_RECURSE CPP_FILES src/*.cpp)
# find *.c
file(GLOB_RECURSE C_FILES src/*.c)
# find *.f90
file(GLOB_RECURSE F90_FILES src/*.f90)
# add all source files
ADD_EXECUTABLE(
axisem3d
${CPP_FILES}
${C_FILES}
${F90_FILES}
# non-header-only local externals
external/yaml/Yaml.cpp
)


################# linking #################
target_link_libraries(
axisem3d
${MPI_C_LIBRARIES}
${FFTW_LIBRARIES}
${METIS_LIBRARIES}
${NETCDF_LIBRARIES}
${HDF5_LIBRARIES}
${HDF5_HL_LIBRARIES}
${ADDITIONAL_LIBS}
)
