#ifndef KINETICENERGYKERNEL_H
#define KINETICENERGYKERNEL_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>


// Host-side function to launch the kernel
void launchKineticEnergyKernel(double3* d_velocities, double* d_masses, double* d_inverseMasses, double3* d_totalForces, double timeStep, int numAtoms, double* d_kineticEnergies, double* d_totalKEnergy);

#endif // KINETICENERGYKERNEL_H