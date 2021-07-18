cmake_minimum_required(VERSION 3.16)
project(cpp_freqtable)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)
message("CMAKE_CXX_EXTENSIONS is " ${CMAKE_CXX_EXTENSIONS})
set(CMAKE_CXX20_STANDARD_COMPILE_OPTION "-std=c++11")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -O3 -march=native")
message( "CMAKE_CXX20_STANDARD_COMPILE_OPTION is " ${CMAKE_CXX20_STANDARD_COMPILE_OPTION} ".")
set(CMAKE_VERBOSE_MAKEFILE on)

add_executable(cpp_freqtable src/main.cpp)
add_executable(readstofreq src/csv_from_reads.cpp)