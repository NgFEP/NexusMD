if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Please run CMake from the root directory.")
endif()

find_package(Python REQUIRED) 
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_SOURCE_DIR}/Nexus/core)
include_directories(${CMAKE_SOURCE_DIR}/Nexus/core/dispatcher/headers)
include_directories(${CMAKE_SOURCE_DIR}/Nexus/core/common/headers)
include_directories(${CMAKE_SOURCE_DIR}/Nexus/core/baseline/headers)
include_directories(${CMAKE_SOURCE_DIR}/Nexus/core/baseline/headers/Forces)
include_directories(${CMAKE_SOURCE_DIR}/Nexus/core/libs/pugixml-1.14/src)

set(CMAKE_SWIG_FLAGS "-py3")

set(SWIG_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Nexus.i)

set_property(SOURCE ${SWIG_INTERFACE} PROPERTY CPLUSPLUS ON)

swig_add_library(Nexus LANGUAGE python SOURCES ${SWIG_INTERFACE})

set_target_properties(Nexus PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
)

add_custom_command(
    TARGET Nexus
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/wrapper"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_BINARY_DIR}/wrapper/Nexus.py"
        "${CMAKE_BINARY_DIR}/wrapper/Release/Nexus.py"
)


#swig_link_libraries(Nexus ${PYTHON_LIBRARIES} baseline common dispatcher libs)
swig_link_libraries(Nexus ${PYTHON_LIBRARIES} pugixml_lib common baseline dispatcher)

add_custom_command(
    TARGET Nexus
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/Nexus/core/libs/fftw-3.3.5-dll64/libfftw3-3.dll"
        "$<TARGET_FILE_DIR:Nexus>"  # This is the folder where Nexus.pyd/.dll goes
)