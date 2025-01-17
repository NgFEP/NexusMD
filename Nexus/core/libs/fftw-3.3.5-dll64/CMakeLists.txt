

add_library(fftw SHARED IMPORTED)
set_target_properties(fftw PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.dll"
    IMPORTED_IMPLIB "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64/libfftw3-3.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/fftw-3.3.5-dll64"
)
