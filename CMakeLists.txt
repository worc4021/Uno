# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
        DESCRIPTION "Uno (Unifying Nonconvex Optimization)"
        LANGUAGES CXX C Fortran)
set(CMAKE_CXX_STANDARD 17)

include(GNUInstallDirs)
set(MKL_LINK static)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL REQUIRED)

include($ENV{CMAKE_INSTALL_PREFIX}/cmake/GKlibTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/metisTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/AMPLTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/bqpdTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/mumpsTargets.cmake)

include(CMakeGraphVizOptions)
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
                            MACRO_NAMESPACE "FC_"
                            SYMBOL_NAMESPACE "FC_")

# directories
set(DIRECTORIES uno ${CMAKE_BINARY_DIR}/include)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")


# source files
file(GLOB_RECURSE UNO_SOURCE_FILES uno/*.cpp)

file(GLOB_RECURSE UNO_HEADER_FILES uno/*.hpp)

set(UNO_HEADER_REL)
foreach(file ${UNO_HEADER_FILES})
   file(RELATIVE_PATH relpath ${CMAKE_CURRENT_SOURCE_DIR} ${file})
   list(APPEND UNO_HEADER_REL ${relpath})
endforeach()

# unit test source files
file(GLOB TESTS_UNO_SOURCE_FILES unotest/*.cpp)

#########################
# external dependencies #
#########################
set(LIBRARIES)
set(UNO_DEFINES)


# HSL or MA57
find_library(HSL hsl)
if(HSL)
   link_to_uno(hsl ${HSL})
else()
   message(WARNING "Optional library HSL was not found.")
   find_library(MA57 ma57)
   if(MA57)
      link_to_uno(ma57 ${MA57})
   else()
      message(WARNING "Optional library MA57 was not found.")
   endif()
endif()
if(HSL OR MA57)
   list(APPEND UNO_SOURCE_FILES uno/solvers/MA57/MA57Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MA57SolverTests.cpp)
   
   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()

list(APPEND LIBRARIES metis::metis)
list(APPEND UNO_DEFINES HAS_METIS)

list(APPEND UNO_SOURCE_FILES uno/solvers/BQPD/BQPDSolver.cpp uno/solvers/BQPD/wdotd.f)
list(APPEND UNO_DEFINES HAS_BQPD)
list(APPEND LIBRARIES bqpd::bqpd_sparse)

list(APPEND UNO_SOURCE_FILES uno/solvers/MUMPS/MUMPSSolver.cpp)
list(APPEND TESTS_UNO_SOURCE_FILES unotest/MUMPSSolverTests.cpp)
list(APPEND LIBRARIES mumps::libdmumps)
list(APPEND UNO_DEFINES MUMPS_SEQUENTIAL)
list(APPEND UNO_DEFINES HAS_MUMPS)

find_package(HIGHS)
if(NOT HIGHS)
   message(WARNING "Optional library HiGHS was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/HiGHS/HiGHSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/HiGHSSolverTests.cpp)
   link_to_uno(highs ${HIGHS})
   list(APPEND LIBRARIES highs::highs)
endif()

# MUMPS
find_package(MUMPS)
if(NOT MUMPS_LIBRARY)
   message(WARNING "Optional library MUMPS was not found.")
else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/MUMPS/MUMPSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MUMPSSolverTests.cpp)
   list(APPEND LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})
   
   list(APPEND DIRECTORIES ${MUMPS_INCLUDE_DIR})


###############
# Uno library #
###############
add_library(uno STATIC ${UNO_SOURCE_FILES})
set_target_properties(uno PROPERTIES 
      POSITION_INDEPENDENT_CODE ON
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_include_directories(uno PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/uno> $<INSTALL_INTERFACE:include/uno>)
target_include_directories(uno PRIVATE ${CMAKE_BINARY_DIR}/include)
target_compile_definitions(uno PUBLIC ${UNO_DEFINES})
target_link_libraries(uno PUBLIC ${LIBRARIES})


target_compile_options(uno PUBLIC -Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
endif()


# copy the option file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_sources(uno PUBLIC
      FILE_SET public_headers
      TYPE HEADERS
      BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
      FILES ${UNO_HEADER_FILES}
)  

######################
# optional AMPL main #
######################

add_executable(uno_ampl bindings/AMPL/AMPLModel.cpp bindings/AMPL/uno_ampl.cpp)

target_link_libraries(uno_ampl PUBLIC uno ampl::ampl)
add_definitions("-D HAS_AMPLSOLVER")
set_target_properties(uno_ampl PROPERTIES 
   POSITION_INDEPENDENT_CODE ON
   MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

##################################
# optional GoogleTest unit tests #
##################################
if(WITH_GTEST)
   find_package(GTest CONFIG REQUIRED)
   if(NOT GTest_DIR)
      message(WARNING "Optional library GTest was not found.")
   else()
      add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
      target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
   endif()
endif()

add_subdirectory(localtest)

# install the library
install(TARGETS uno
      EXPORT unoTargets
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${CMAKE_BUILD_TYPE}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${CMAKE_BUILD_TYPE}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${CMAKE_BUILD_TYPE}
      FILE_SET public_headers)

install(EXPORT unoTargets
      FILE unoTargets.cmake
      NAMESPACE uno::
      DESTINATION cmake)