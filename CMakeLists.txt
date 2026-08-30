cmake_minimum_required(VERSION 3.16)
project(lattice LANGUAGES Fortran)

add_library(
    lattice
        src/lattice_mod.f90
)

target_compile_options(
    lattice
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion-extra
        -Wimplicit-interface
        -fcheck=all
        -fbacktrace
        -ffpe-trap=invalid,zero,overflow
        -g
        -O0
        -ffree-line-length-none
)

add_executable(
    main
        examples/main.f90   
)

target_link_libraries(
    main
    PRIVATE
        lattice
)