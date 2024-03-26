# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
			DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
            LANGUAGES CXX C Fortran)

# set C++17 and enable other languages
set(CMAKE_CXX_STANDARD 20)
enable_language(CXX C Fortran)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion -Wmaybe-uninitialized") # -fsanitize=address
    set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion -Wmaybe-uninitialized")
    set(CMAKE_CXX_FLAGS_DEBUG "-pg")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if(WIN32)
        set(CMAKE_CXX_FLAGS "/Wall /Wextra /Wnon-virtual-dtor /Wunused-value /Wunused-function /Wconversion /Wmaybe-uninitialized /fp:precise /EHa")
    else()
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -Wunused-value -Wunused-function -Wconversion -Wmaybe-uninitialized -fp-model=precise")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion -Wmaybe-uninitialized /EHa")
endif()
add_compile_definitions($<$<NOT:$<CONFIG:Debug>>:NDEBUG>)


# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")

option(HAVE_BQPD "Enable BQPD" ON)
option(HAVE_MA27 "Enable MA27" ON)
option(HAVE_MA57 "Enable MA57" OFF)
option(HAVE_AMPL "Enable AMPL" OFF)
option(USE_OPENMP "Enable OpenMP" OFF)

# directories
set(DIRECTORIES uno)

# source files
file(GLOB_RECURSE UNO_SOURCE_FILES
    uno/Uno.cpp
    uno/ingredients/*.cpp
    uno/optimization/*.cpp
    uno/preprocessing/*.cpp
    uno/tools/*.cpp
)

# find libraries
set(LIBRARIES)
set(DEFINES)
set(OPTIONAL_LIBRARIES amplsolver ma57 metis bqpd CACHE STRING "Optional libraries")
set(REQUIRED_LIBRARIES dl blas lapack)

add_subdirectory(resources)

if (HAVE_BQPD)
list(APPEND UNO_SOURCE_FILES uno/solvers/QP/BQPDSolver.cpp uno/solvers/QP/wdotd.f)
list(APPEND LIBRARIES bqpd::bqpd)
list(APPEND DEFINES HAS_BQPD)
endif()

if (HAVE_MA27)
list(APPEND UNO_SOURCE_FILES uno/solvers/linear/MA27Solver.cpp)
list(APPEND LIBRARIES hsl::ma27)
list(APPEND DEFINES HAS_MA27)
endif()

if (HAVE_MA57)
list(APPEND UNO_SOURCE_FILES uno/solvers/linear/MA57Solver.cpp)
list(APPEND LIBRARIES hsl::ma57)
list(APPEND DEFINES HAS_MA57)
endif()

if (HAVE_AMPL)
list(APPEND UNO_SOURCE_FILES uno/interfaces/AMPL/AMPLModel.cpp)
list(APPEND LIBRARIES amplsolver)
list(APPEND DEFINES HAS_AMPL)
endif()

if (USE_OPENMP)
find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)
endif()

#######
# Uno #
#######
add_library(uno STATIC ${UNO_SOURCE_FILES})
add_library(Uno::uno ALIAS uno)
target_include_directories(uno PUBLIC ${DIRECTORIES})
target_compile_definitions(uno PUBLIC ${DEFINES})

set_property(TARGET uno PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# link the libraries
target_link_libraries(uno PUBLIC ${LIBRARIES})


# copy the option file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#############
# AMPL main #
#############
if(HAVE_AMPL)
add_executable(uno_ampl uno/main.cpp)
target_link_libraries(uno_ampl PUBLIC uno)
endif()

#########################
# GoogleTest unit tests #
#########################
if(WITH_GTEST)
    find_package(GTest CONFIG REQUIRED)
    if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
        file(GLOB TESTS_UNO_SOURCE_FILES
            unotest/*.cpp
        )
        add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
        target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
    endif()
endif()

add_subdirectory(nlp_test)