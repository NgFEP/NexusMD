# CMake file for the CUDA module
cmake_minimum_required(VERSION 3.18)

# Define the CUDA module
project(cuda LANGUAGES CXX CUDA)

# Enable relocatable device code (RDC)
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -rdc=true")

# Enable relocatable device code (RDC)
#set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -rdc=true")

# Define supported CUDA architectures
# set(CUDA_ARCHITECTURES "60;70;75;80;86;89")

# Source and header files
file(GLOB CUDA_HEADERS 
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces/kernels/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/kernels/*.h"
)
file(GLOB CUDA_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Forces/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Forces/kernels/*.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/kernels/*.cu"
)

# Add the library
add_library(cuda STATIC ${CUDA_HEADERS} ${CUDA_SOURCES})

# Include directories
target_include_directories(cuda PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/headers"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/Forces/kernels"
    "${CMAKE_CURRENT_SOURCE_DIR}/headers/kernels"
    "${CMAKE_SOURCE_DIR}/core/common/headers"
    "${CMAKE_SOURCE_DIR}/core/libs/pugixml-1.14/src"    # Pugixml headers
    "${CMAKE_SOURCE_DIR}/core/libs/fftw-3.3.5-dll64"    # FFTW headers
)


if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    message(FATAL_ERROR "CMAKE_CUDA_ARCHITECTURES is not defined in the cuda module.")
endif()
message(STATUS "CUDA Subdirectory: CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}")

# Optional: Print diagnostic information
message(STATUS "CUDA Architectures in CUDA module: ${CMAKE_CUDA_ARCHITECTURES}")

set_target_properties(cuda PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

# Compiler options for CUDA
target_compile_options(cuda PRIVATE 
    $<$<COMPILE_LANGUAGE:CUDA>:--generate-line-info -use_fast_math -rdc=true>
)

# Find and link CUDA Toolkit
find_package(CUDAToolkit REQUIRED)
find_library(CUDA_CUDADEVRT_LIBRARY
    NAMES cudadevrt.lib
    HINTS "${CUDAToolkit_LIBRARY_DIR}"
)

target_link_libraries(cuda PRIVATE
    CUDA::cudart              # CUDA runtime
    CUDA::nvrtc               # NVRTC
    CUDA::cuda_driver         # CUDA Driver API
    ${CUDA_CUDADEVRT_LIBRARY} # Device runtime (last in order)
    common                    # Common module
    libs                      # FFTW, pugixml, etc.

)

# Enable relocatable device code for the CUDA target


# Compiler options for CUDA
#target_compile_options(cuda PRIVATE 
#    $<$<COMPILE_LANGUAGE:CUDA>:--generate-line-info -use_fast_math -maxrregcount=0>
#)


#target_compile_options(cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-rdc=true>)


#target_link_libraries(cuda PRIVATE ${CUDAToolkit_LIBRARY_DIR}/cudadevrt.lib)

# Optional: Print diagnostic information
message(STATUS "CUDA Architectures: ${CUDA_ARCHITECTURES}")
