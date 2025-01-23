#ifndef VERLET_INTEGRATION_KERNEL_H
#define VERLET_INTEGRATION_KERNEL_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

// Function to launch the Verlet integration kernel
void launchVerletIntegrationKernel(double3* d_atomPositions, double3* d_velocities, double3* d_forces, double* d_inverseMasses, double3* d_boxSize, double3* d_lb, double3* d_ub, int numAtoms, double dt);


#endif // VERLET_INTEGRATION_KERNEL_H
