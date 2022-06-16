cmake_minimum_required(VERSION 3.14)
project(DISTOPIA LANGUAGES C CXX)
enable_testing()

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)


add_compile_options(-O3 -mtune=native -march=native)


add_executable(d_array d_array.cc)
add_executable(batched batched.cc)

# Externally provided libraries
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.10.x)

FetchContent_Declare(googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG master) # need master for benchmark::benchmark

FetchContent_MakeAvailable(
        googletest
        googlebenchmark)

add_executable(bench bench.cc)
target_link_libraries(bench benchmark::benchmark)
