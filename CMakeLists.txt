cmake_minimum_required(VERSION 3.1)
project(pde_solver)

include(CheckCXXCompilerFlag)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion -Wsign-conversion")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion")
    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if (HAS_CPP14_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    else()
        message(FATAL_ERROR "Unsupported compiler -- pde_solver requires C++14 support!")
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

set(PDE_SOLVER_SRCS
    closed_form.hpp
    closed_form.cpp
    solver.hpp
    solver.cpp
    payoff.hpp
    payoff.cpp
    volatility.hpp
    volatility.cpp
    boundary_conditions.hpp
    boundary_conditions.cpp
    rate.hpp 
    rate.cpp
    mesh.hpp
    mesh.cpp
    matrix.hpp
    matrix.cpp
    pde.hpp
    pde.cpp
    greeks.hpp
    greeks.cpp
    main.cpp
)

set(PDE_SOLVER_TARGET pde_solver)
add_executable(${PDE_SOLVER_TARGET} ${PDE_SOLVER_SRCS})

