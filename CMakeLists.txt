cmake_minimum_required(VERSION 3.18)

project(NexusMD LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_SWIG_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper)


add_subdirectory(NexusMD)
add_subdirectory(wrapper)


include_directories(
    ${CMAKE_SOURCE_DIR}/NexusMD/core/libs/fftw-3.3.5-dll64
    ${CMAKE_SOURCE_DIR}/NexusMD/core/libs/pugixml-1.14/src
)

link_directories(
    ${CMAKE_SOURCE_DIR}/NexusMD/core/libs/pugixml-1.14/build/Release
    ${CMAKE_SOURCE_DIR}/NexusMD/core/libs/fftw-3.3.5-dll64
)


target_link_options(NexusMD PRIVATE /INCREMENTAL)

set_target_properties(NexusMD PROPERTIES LINK_LIBRARY_DEPENDENCIES ON)

target_link_libraries(NexusMD FFTW3::fftw3 pugixml_lib)
