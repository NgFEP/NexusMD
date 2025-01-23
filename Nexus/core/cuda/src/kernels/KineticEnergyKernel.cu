#include "KineticEnergyKernel.h"


// GPU kernel to calculate kinetic energy for each atom
__global__ void KineticEnergyKernel(double3* d_velocities, double* d_masses, double* d_inverseMasses, double3* d_totalForces, double timeShift, int numAtoms, double* d_kineticEnergies, double* d_totalKEnergy) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= numAtoms) return;  // Ensure we don't access out-of-numAtoms
    double mass = d_masses[idx];
    double invMass = d_inverseMasses[idx];

    double3 velocity = d_velocities[idx];
    double3 force = d_totalForces[idx];

    // Apply time shift to velocity
    velocity.x += force.x * invMass * timeShift;
    velocity.y += force.y * invMass * timeShift;
    velocity.z += force.z * invMass * timeShift;

    // Compute kinetic energy with shifted velocity
    d_kineticEnergies[idx] = 0.5 * mass * (velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z);

    // Write the result to the kinetic energy array
    //d_totalKEnergy += d_kineticEnergies[idx];
    atomicAdd(d_totalKEnergy, d_kineticEnergies[idx]);



}

// Host function to launch the KineticEnergyKernel
void launchKineticEnergyKernel(double3* d_velocities, double* d_masses, double* d_inverseMasses, double3* d_totalForces, double timeStep, int numAtoms, double* d_kineticEnergies, double* d_totalKEnergy) {
    // Calculate half time step for time shift
    double timeShift = 0.5 * timeStep;


    // Define block size and number of blocks
    int blockSize = 256;
    int numBlocks = (numAtoms + blockSize - 1) / blockSize;



    // Launch the CUDA kernel
    KineticEnergyKernel <<<numBlocks, blockSize >>> (d_velocities, d_masses, d_inverseMasses, d_totalForces, timeShift, numAtoms, d_kineticEnergies, d_totalKEnergy);

    // Synchronize the device
    //cudaDeviceSynchronize();




}
