cmake_minimum_required(VERSION 3.18)

project(baseline)

file(GLOB BASELINE_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces/*.h"
)

file(GLOB BASELINE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Forces/*.cpp"
)

add_library(baseline STATIC ${BASELINE_HEADERS} ${BASELINE_SOURCES})

target_include_directories(baseline PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/headers"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces"
    "${CMAKE_SOURCE_DIR}/Nexus/core/common/headers"
    "${CMAKE_SOURCE_DIR}/Nexus/core/libs/pugixml-1.14/src"
)

target_link_libraries(baseline PRIVATE 
    pugixml_lib
    FFTW3::fftw3
    common
)
