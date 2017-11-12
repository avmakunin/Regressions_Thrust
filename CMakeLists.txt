cmake_minimum_required(VERSION 3.8)
project(Regressions_Thrust)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES src/main.cu src/koeffs.h)
add_executable(Regressions_Thrust ${SOURCE_FILES})