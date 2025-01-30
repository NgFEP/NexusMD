cmake_minimum_required(VERSION 3.18)

project(common)

file(GLOB COMMON_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h"
)
file(GLOB COMMON_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_library(common STATIC ${COMMON_HEADERS} ${COMMON_SOURCES})

target_include_directories(common PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/headers"
)

find_package(CUDAToolkit REQUIRED)

target_link_libraries(common PRIVATE 
    CUDA::cudart 
    pugixml_lib 
)

set_target_properties(common PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_STANDARD 11
    CUDA_STANDARD_REQUIRED ON
)
