cmake_minimum_required(VERSION 3.26)
project(MANGLING C Fortran)

add_executable(xintface Fintface.f Cintface.c)
