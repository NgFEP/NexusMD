#include "VerletIntegrationKernel.h"


//verlet inegration

// Device function to apply periodic boundary conditions (PBC)
__device__ void applyPBC(
    double3* d_atomPositions,
    double3* d_boxSize,
    double3* d_lb,
    double3* d_ub) 
{
    if (d_atomPositions->x < d_lb->x) {
        d_atomPositions->x += d_boxSize->x;
    }
    else if (d_atomPositions->x > d_ub->x) {
        d_atomPositions->x -= d_boxSize->x;
    }

    if (d_atomPositions->y < d_lb->y) {
        d_atomPositions->y += d_boxSize->y;
    }
    else if (d_atomPositions->y > d_ub->y) {
        d_atomPositions->y -= d_boxSize->y;
    }

    if (d_atomPositions->z < d_lb->z) {
        d_atomPositions->z += d_boxSize->z;
    }
    else if (d_atomPositions->z > d_ub->z) {
        d_atomPositions->z -= d_boxSize->z;
    }
}


// Kernel to update positions using Verlet integration and apply PBC
__global__ void verletIntegrationKernel(
    double3* d_atomPositions,
    double3* d_velocities,
    double3* d_forces,
    double* d_inverseMasses,
    double3* d_boxSize,
    double3* d_lb,
    double3* d_ub,
    int numAtoms,
    double dt)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= numAtoms) return;  // Ensure we don't access out-of-bounds

    if (d_inverseMasses[idx] != 0.0) {
        // Update velocity (half-step)
        d_velocities[idx].x += d_inverseMasses[idx] * d_forces[idx].x * dt;
        d_velocities[idx].y += d_inverseMasses[idx] * d_forces[idx].y * dt;
        d_velocities[idx].z += d_inverseMasses[idx] * d_forces[idx].z * dt;

        // Update position (full-step)
        d_atomPositions[idx].x += d_velocities[idx].x * dt;
        d_atomPositions[idx].y += d_velocities[idx].y * dt;
        d_atomPositions[idx].z += d_velocities[idx].z * dt;

        // Apply periodic boundary conditions (PBC)
        applyPBC(d_atomPositions, d_boxSize, d_lb, d_ub);


        // Store updated positions back to global memory
    }
}

void launchVerletIntegrationKernel(
    double3* d_atomPositions,
    double3* d_velocities,
    double3* d_forces,
    double* d_inverseMasses,
    double3* d_boxSize,
    double3* d_lb,
    double3* d_ub,
    int numAtoms,
    double dt
) {

    // Define block size and number of blocks
    int blockSize = 256;
    int numBlocks = (numAtoms + blockSize - 1) / blockSize;

    // Launch the kernel
    verletIntegrationKernel <<<numBlocks, blockSize >>> (d_atomPositions, d_velocities, d_forces, d_inverseMasses,d_boxSize, d_lb, d_ub, numAtoms, dt);
}



