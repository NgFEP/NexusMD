#include "InverseMassesKernel.h"

//calculating d_inverseMasses

// Kernel function to compute inverse masses
__global__ void computeInverseMassesKernel(double* d_masses, double* d_inverseMasses, int numAtoms) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Check if the thread index is within bounds
    if (idx >= numAtoms) return;

    // Compute the inverse mass
    if (d_masses[idx] == 0.0) {
        d_inverseMasses[idx] = 0.0;
    }
    else {
        d_inverseMasses[idx] = 1.0 / d_masses[idx];
    }
}

// Function to launch the kernel
void launchComputeInverseMassesKernel(double* d_masses, double* d_inverseMasses, int numAtoms) {
    int blockSize = 256;  // Number of threads per block
    int numBlocks = (numAtoms + blockSize - 1) / blockSize;  // Calculate the number of blocks

    // Launch the kernel
    computeInverseMassesKernel <<<numBlocks, blockSize >>> (d_masses, d_inverseMasses, numAtoms);

    // Check for errors in kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA kernel launch failed: %s\n", cudaGetErrorString(err));
    }
}


