# CMake file for the libs project
cmake_minimum_required(VERSION 3.18)

project(libs)  # Create a project for libs

# Add an INTERFACE library for libs
add_library(libs INTERFACE)

# Add FFTW
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64")
    target_include_directories(libs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64")
    target_link_libraries(libs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.lib")
else()
    message(WARNING "FFTW directory not found: ${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64")
endif()

# Add pugixml
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/build/Release/pugixml.lib")
    target_include_directories(libs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/src")
    target_link_libraries(libs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/build/Release/pugixml.lib")
else()
    message(WARNING "Pugixml directory not found: ${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14")
endif()
