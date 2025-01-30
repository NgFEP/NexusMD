cmake_minimum_required(VERSION 3.18)
project(libs)

# ------------------------------------------------------------
# 1) PugiXML
# ------------------------------------------------------------
add_library(pugixml_lib INTERFACE)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/build/Release/pugixml.lib")
    target_include_directories(pugixml_lib INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/src")
    target_link_libraries(pugixml_lib INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14/build/Release/pugixml.lib")
else()
    message(WARNING "Pugixml library not found in: ${CMAKE_CURRENT_SOURCE_DIR}/pugixml-1.14")
endif()

# ------------------------------------------------------------
# 2) FFTW3
# ------------------------------------------------------------

add_library(FFTW3::fftw3 SHARED IMPORTED GLOBAL)


set_target_properties(FFTW3::fftw3 PROPERTIES
    IMPORTED_IMPLIB "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.lib"
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.dll"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64"
)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.lib")
    message(WARNING "FFTW3 library not found in: ${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64")
endif()
