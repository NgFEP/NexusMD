cmake_minimum_required(VERSION 3.18)

# Enforce running from top-level
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Please run CMake from the root directory.")
endif()

find_package(Python REQUIRED) 
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(
    "${CMAKE_SOURCE_DIR}/Nexus/core"
    "${CMAKE_SOURCE_DIR}/Nexus/core/dispatcher/headers"
    "${CMAKE_SOURCE_DIR}/Nexus/core/common/headers"
    "${CMAKE_SOURCE_DIR}/Nexus/core/libs/pugixml-1.14/src"
    "${CMAKE_SOURCE_DIR}/Nexus/core/libs/fftw-3.3.5-dll64"
)

# Set SWIG flags
set(CMAKE_SWIG_FLAGS "-py3")

# Define SWIG interface file
set(SWIG_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Nexus.i)
set_property(SOURCE ${SWIG_INTERFACE} PROPERTY CPLUSPLUS ON)

set(DUMMY_CU "${CMAKE_CURRENT_SOURCE_DIR}/dummy.cu")
# Mark the dummy file as CUDA
set_source_files_properties("${DUMMY_CU}" PROPERTIES LANGUAGE CUDA)

# Add the SWIG-generated library
swig_add_library(Nexus LANGUAGE python SOURCES "${SWIG_INTERFACE}" "${DUMMY_CU}")

# Set CUDA properties for the wrapper target
set_target_properties(Nexus PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/wrapper
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
    LINKER_LANGUAGE CUDA
)

target_compile_options(Nexus PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>
)

find_package(CUDAToolkit REQUIRED)

# Locate cudadevrt
#find_library(CUDA_CUDADEVRT_LIBRARY
#    NAMES cudadevrt.lib
#    HINTS "${CUDAToolkit_LIBRARY_DIR}"
#    PATH_SUFFIXES "stubs"
#)
if(NOT CUDA_CUDADEVRT_LIBRARY)
    message(FATAL_ERROR "Could not find cudadevrt library.")
endif()

swig_link_libraries(Nexus
    cuda                # Custom CUDA library
    dispatcher          # Dispatcher module
    common              # Common module
    baseline            # Baseline module
    ${PYTHON_LIBRARIES} # Python library
    #${CUDA_CUDADEVRT_LIBRARY} # cudadevrt last
)

# Ensure the `cuda` target is fully built before the `wrapper` target
add_dependencies(Nexus cuda)
