cmake_minimum_required(VERSION 3.18)

project(dispatcher)

file(GLOB DISPATCHER_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h"
)
file(GLOB DISPATCHER_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_library(dispatcher STATIC ${DISPATCHER_HEADERS} ${DISPATCHER_SOURCES})

target_include_directories(dispatcher PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/headers"
    "${CMAKE_SOURCE_DIR}/NexusMD/core/common/headers"


)

target_link_libraries(dispatcher PRIVATE 
    pugixml_lib
    common        
    baseline
    
)

