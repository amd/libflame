##Copyright (C) 2020, Advanced Micro Devices, Inc.## 

enable_language(Fortran)

option (BUILD_SINGLE "enabled single testing variant for netlib" ON)
option (BUILD_COMPLEX "enabled complex testing variant for netlib" ON)
option (BUILD_DOUBLE "enabled double testing variant for netlib" ON)
option (BUILD_COMPLEX16 "enabled complex16 testing variant for netlib" ON)

add_subdirectory(netlib_win_test)
