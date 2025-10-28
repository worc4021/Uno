cmake_minimum_required(VERSION 3.25)
# define the project
project(Uno VERSION 2.3.0
        DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
        LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
                            MACRO_NAMESPACE "FC_"
                            SYMBOL_NAMESPACE "FC_")
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")
option(WITH_MUMPS "Enable MUMPS support" ON)
message(STATUS "MUMPS Solver: WITH_MUMPS=${WITH_MUMPS}")
option(WITH_BQPD "Enable BQPD support" ON)
message(STATUS "BQPD Solver: WITH_BQPD=${WITH_BQPD}")
option(WITH_AMPL "Enable AMPL support" ON)
message(STATUS "AMPL Interface: WITH_AMPL=${WITH_AMPL}")

include(FetchContent)

if (WITH_GTEST)
enable_testing()
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG v1.17.0
        EXCLUDE_FROM_ALL
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
endif(WITH_GTEST)

# source files
file(GLOB UNO_SOURCE_FILES
   interfaces/C/Uno_C_API.cpp
   uno/Uno.cpp
   uno/ingredients/constraint_relaxation_strategies/*.cpp
   uno/ingredients/globalization_mechanisms/*.cpp
   uno/ingredients/globalization_strategies/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/filters/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/funnel_methods/*.cpp
   uno/ingredients/hessian_models/*.cpp
   uno/ingredients/inequality_handling_methods/*.cpp
   uno/ingredients/inequality_handling_methods/inequality_constrained_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/*.cpp
   uno/ingredients/inertia_correction_strategies/*.cpp
   uno/ingredients/subproblem/*.cpp
   uno/ingredients/subproblem_solvers/*.cpp
   uno/model/*.cpp
   uno/optimization/*.cpp
   uno/options/*.cpp
   uno/tools/*.cpp
)
add_library(uno ${UNO_SOURCE_FILES})
target_include_directories(uno PUBLIC uno)
target_include_directories(uno PRIVATE ${CMAKE_BINARY_DIR}/include)

if (WITH_MUMPS)
        target_compile_definitions(uno PRIVATE -DWITH_MUMPS)
        
        FetchContent_Declare(
                mumps
                GIT_REPOSITORY https://github.com/worc4021/MUMPS_cmake.git
                GIT_TAG main
                EXCLUDE_FROM_ALL
        )
        FetchContent_MakeAvailable(mumps)

        get_target_property(UNO_SOURCE_FILES uno SOURCES)
        target_link_libraries(uno PUBLIC libdmumps)

        list(APPEND UNO_SOURCE_FILES
                ${CMAKE_CURRENT_SOURCE_DIR}/uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp
        )
        
        set_target_properties(uno PROPERTIES SOURCES "${UNO_SOURCE_FILES}")
endif()

if (WITH_BQPD)
    target_compile_definitions(uno PRIVATE -DWITH_BQPD)
    set(BQPD_INTERFACE OFF)
    FetchContent_Declare(
                bqpd
                GIT_REPOSITORY git@github.com:worc4021/upgraded-octo-adventure.git
                GIT_TAG master
                EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(bqpd)

    get_target_property(UNO_SOURCE_FILES uno SOURCES)
    list(APPEND UNO_SOURCE_FILES
           uno/ingredients/subproblem_solvers/BQPD/BQPDSolver.cpp
           uno/ingredients/subproblem_solvers/BQPD/BQPDEvaluationSpace.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/BQPDSolverTests.cpp)
   set_target_properties(uno PROPERTIES SOURCES "${UNO_SOURCE_FILES}")
   target_link_libraries(uno PUBLIC bqpd_sparse)
endif()

if (WITH_AMPL)
    target_compile_definitions(uno PRIVATE -DWITH_AMPL)
    
    FetchContent_Declare(
        amplasl
        GIT_REPOSITORY https://github.com/worc4021/AMPL.git
        GIT_TAG main
        EXCLUDE_FROM_ALL
        )
FetchContent_MakeAvailable(amplasl)
add_executable(uno_ampl interfaces/AMPL/AMPLModel.cpp interfaces/AMPL/uno_ampl.cpp)
target_link_libraries(uno_ampl PRIVATE uno ampl)
endif()