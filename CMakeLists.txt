cmake_minimum_required(VERSION 3.18)
project(Nexus LANGUAGES CXX CUDA)

# Enable CUDA and C++
enable_language(CXX)
enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Optionally define CUDA architectures globally
# so all subprojects share them:
#set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86;89")
#set(CMAKE_CUDA_ARCHITECTURES "89")
#set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86;89" CACHE STRING "CUDA architectures")
#set(CMAKE_CUDA_ARCHITECTURES "89" CACHE STRING "CUDA architectures")
#set(CMAKE_CUDA_ARCHITECTURES "89" CACHE STRING "CUDA architectures for all subprojects")
set(CMAKE_CUDA_ARCHITECTURES "89" CACHE STRING "CUDA architectures for all subprojects" FORCE)
message(STATUS "Root CMakeLists: CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}")

# Optionally also force dynamic runtime:
set(CMAKE_CUDA_RUNTIME_LIBRARY "Shared" CACHE STRING "" FORCE)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_SWIG_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper)

#set(CUDA_ARCHITECTURES_PROPAGATED "${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "Root CMakeLists: CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}")

add_subdirectory(Nexus)
add_subdirectory(wrapper)
