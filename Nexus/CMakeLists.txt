# Add the core subdirectory
add_subdirectory(core)

# Propagate CUDA_ARCHITECTURES explicitly to subdirectories
set(CMAKE_CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}" PARENT_SCOPE)
message(STATUS "Nexus Subdirectory: CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}")
