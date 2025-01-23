#ifndef INVERSE_MASSES_KERNEL_H
#define INVERSE_MASSES_KERNEL_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

void launchComputeInverseMassesKernel(double* d_masses, double* d_inverseMasses, int numAtoms);


#endif // INVERSE_MASSES_KERNEL_H
