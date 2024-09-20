#include <cuda.h>
#include <cuda_runtime.h> // Stops underlining of __global__


// Declaration of the device function
__device__ void minimumImageVectorDevice(double3* coords1, double3* coords2, double3* delta, double3* d_boxsize);
