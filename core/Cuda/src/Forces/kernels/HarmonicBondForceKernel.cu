// HarmonicBondForce.cu
#include <PeriodicBoundaryConditionKernel.h>
#include "HarmonicBondForceKernel.h"
#include <vector>

//. using a  __device__ const int blocksize = 320; Is faster than using blockDim.x all the time
//__device__ const int blockSize = 256;

//__global__ void BondForcesKernel_shared(
//    const Coords3D* atomPositions,
//    const HBondParams* bondParams,
//    Coords3D* forces,
//    double* totalPEnergy,
//    int numBonds,
//    const PeriodicBoundaryCondition::BoxInfo boxInfo)
//{
//    // Allocate shared memory for a block of atom positions and bond parameters
//
//
//
//
//    extern __shared__ Coords3D sharedAtomPositions[];
//    __shared__ HBondParams sharedBondParams[blockSize];  // Assuming max 256 threads/block
//    __shared__ double sharedEnergy[1];
//
//    int idx = threadIdx.x + blockIdx.x * blockDim.x;
//
//    if (idx >= numBonds) return;  // Ensure we don't access out-of-bounds
//
//    // Load bond parameters into shared memory
//    sharedBondParams[threadIdx.x] = bondParams[idx];
//    __syncthreads();
//
//    // Load atom positions into shared memory
//    const Coords3D& pos1 = atomPositions[sharedBondParams[threadIdx.x].p1];
//    const Coords3D& pos2 = atomPositions[sharedBondParams[threadIdx.x].p2];
//    sharedAtomPositions[threadIdx.x * 2] = pos1;
//    sharedAtomPositions[threadIdx.x * 2 + 1] = pos2;
//    __syncthreads();
//
//    // Load bond parameters for the current bond from shared memory
//    const HBondParams& params = sharedBondParams[threadIdx.x];
//    const Coords3D& sPos1 = sharedAtomPositions[threadIdx.x * 2];
//    const Coords3D& sPos2 = sharedAtomPositions[threadIdx.x * 2 + 1];
//
//    // Compute the vector from pos1 to pos2
//    Coords3D delta = PeriodicBoundaryCondition::minimumImageVector(sPos2, sPos1, boxInfo);
//
//    // Calculate the distance between the two particles
//    double r = sqrt(delta.dot(delta));
//
//    // Compute the energy contribution
//    double deltaIdeal = r - params.d;
//    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
//
//    // Accumulate energy to totalPEnergy using atomicAdd
//    atomicAdd(totalPEnergy, energy);
//
//    // Compute the derivative of the energy with respect to the distance
//    double dEdR = params.k * deltaIdeal;
//
//    // Normalize the delta vector and scale by dEdR
//    Coords3D force;
//    if (r > 0) {
//        force = delta * (dEdR / r);
//    }
//    else {
//        force = Coords3D(0, 0, 0);  // Prevent division by zero
//    }
//
//    // Update forces using atomicAdd to avoid race conditions
//    atomicAdd(&forces[params.p1], force);
//
//    atomicAdd(&forces[params.p2], -force);
//}


//__device__ double atomicAdd_double(double* address, double val) {
//    unsigned long long int* address_as_ull = (unsigned long long int*)address;
//    unsigned long long int old = *address_as_ull, assumed;
//
//    do {
//        assumed = old;
//        old = atomicCAS(address_as_ull, assumed,
//            __double_as_longlong(val + __longlong_as_double(assumed)));
//    } while (assumed != old);
//
//    return __longlong_as_double(old);
//}



__global__ void BondForcesKernel_global(
    double3* d_atomPositions,
    HBondParams* d_bondParams,
    double3* d_forces,
    double* d_totalPEnergy,
    int d_numBonds,
    double3* d_boxsize)  // Periodic boundary conditions info
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= d_numBonds) return;  // Ensure we don't access out-of-bounds

    // Load bond parameters from global memory
    HBondParams& params = d_bondParams[idx];

    // Load atom positions from global memory
    double3 coords1 = d_atomPositions[params.p1];
    double3 coords2 = d_atomPositions[params.p2];
    //printf("params.p1: %f, params.p2: %f\n", params.p2, params.p2);
    //printf("coords1.x: %f, coords1.y: %f, coords1.z: %f\n", coords1.x, coords1.y, coords1.z);

    // Compute the vector from pos1 to pos2 considering periodic boundary conditions
    double3 delta;
    minimumImageVectorDevice(&coords1, &coords2, &delta, d_boxsize);

    // Calculate the distance between the two particles
    double r = sqrt(delta.x*delta.x + delta.y * delta.y + delta.z * delta.z);

    // Compute the energy contribution
    double deltaIdeal = r - params.d;
    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;

    // Accumulate energy to totalPEnergy using atomicAdd
    atomicAdd(d_totalPEnergy, energy);  // You can switch back to atomicAdd if needed

    // Compute the derivative of the energy with respect to the distance
    double dEdR = params.k * deltaIdeal;

    // Normalize the delta vector and scale by dEdR
    double3 force;
    if (r > 0) {
        force.x = delta.x * (dEdR / r);
        force.y = delta.y * (dEdR / r);
        force.z = delta.z * (dEdR / r);
    }
    else {
        force = make_double3(0.0, 0.0, 0.0);  // Prevent division by zero
    }

    // Update forces for both atoms using atomicAdd to avoid race conditions
    //atomicAdd(&forces[params.p1], force);
    //atomicAdd(&forces[params.p2], -force);

    //atomicAddCoords(forces[params.p1], force);
    //atomicAddCoords(forces[params.p2], Coords3D{ -force[0], -force[1], -force[2]});
    // Update forces for both atoms using atomicAdd to avoid race conditions


    atomicAdd(&d_forces[params.p1].x, force.x);
    atomicAdd(&d_forces[params.p1].y, force.y);
    atomicAdd(&d_forces[params.p1].z, force.z);
    atomicAdd(&d_forces[params.p2].x, -force.x);
    atomicAdd(&d_forces[params.p2].y, -force.y);
    atomicAdd(&d_forces[params.p2].z, -force.z);
}


//correct for double
void launchKernelBondForces(
    double3* atomPositions,  // Host-side Coords3D
    HBondParams* bondParams,
    double3* forces,
    double* totalPEnergy,
    double3* boxInfo,
    int numAtoms,
    int numBonds
) {

    // Device memory allocation
    double3* d_atomPositions;
    HBondParams* d_bondParams;
    double3* d_forces;
    double* d_totalPEnergy;
    double3* d_boxsize;


    cudaMalloc(&d_atomPositions, numAtoms * sizeof(double3));
    cudaMalloc(&d_bondParams, numBonds * sizeof(HBondParams));
    cudaMalloc(&d_forces, numAtoms * sizeof(double3));
    cudaMalloc(&d_totalPEnergy, sizeof(double));
    cudaMalloc(&d_boxsize, sizeof(double3));

    // Copy data from host to device
    cudaMemcpy(d_atomPositions, atomPositions, numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_bondParams, bondParams, numBonds * sizeof(HBondParams), cudaMemcpyHostToDevice);
    cudaMemcpy(d_forces, forces, numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_totalPEnergy, totalPEnergy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_boxsize, boxInfo, sizeof(double3), cudaMemcpyHostToDevice);

    // Define block size and number of blocks
    int blockSize = 256; // Example block size
    int numBlocks = (numBonds + blockSize - 1) / blockSize;

    // Launch the kernel using double3 types on the device
    BondForcesKernel_global <<< numBlocks, blockSize >>> (d_atomPositions, d_bondParams, d_forces, d_totalPEnergy, numBonds, d_boxsize);
    //BondForcesKernel_global2 << < numBlocks, blockSize >> > (d_atomPositions, d_bondParams,d_forces, d_totalPEnergy, numBonds, d_boxsize);

    // Copy results back to host
    cudaMemcpy(forces, d_forces, numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
    cudaMemcpy(totalPEnergy, d_totalPEnergy, sizeof(double), cudaMemcpyDeviceToHost);


    // Free device memory
    cudaFree(d_atomPositions);
    cudaFree(d_bondParams);
    cudaFree(d_forces);
    cudaFree(d_totalPEnergy);
    cudaFree(d_boxsize);
}









//__global__ void BondForcesKernel_global2(
//    Coords3D* d_atomPositions,
//    HBondParams* d_bondParams,
//    Coords3D* d_forces,
//    double* d_totalPEnergy,
//    int d_numBonds,
//    Coords3D* d_boxsize)  // Periodic boundary conditions info
//{
//    int idx = threadIdx.x + blockIdx.x * blockDim.x;
//
//    if (idx >= d_numBonds) return;  // Ensure we don't access out-of-bounds
//
//    // Load bond parameters from global memory
//    HBondParams& params = d_bondParams[idx];
//
//    // Load atom positions from global memory
//    Coords3D coords1 = d_atomPositions[params.p1];
//    Coords3D coords2 = d_atomPositions[params.p2];
//    //printf("params.p1: %f, params.p2: %f\n", params.p2, params.p2);
//    //printf("coords1.x: %f, coords1.y: %f, coords1.z: %f\n", coords1.x, coords1.y, coords1.z);
//
//    // Compute the vector from pos1 to pos2 considering periodic boundary conditions
//    Coords3D delta;
//
//    delta[0] = coords2[0] - coords1[0];
//    delta[1] = coords2[1] - coords1[1];
//    delta[2] = coords2[2] - coords1[2];
//    Coords3D dd_boxsize = *d_boxsize;
//    if (delta[0] > 0.5 * (dd_boxsize[0])) {
//        delta[0] -= dd_boxsize[0];
//    }
//    else if (delta[0] < -0.5 * dd_boxsize[0]) {
//        delta[0] += dd_boxsize[0];
//    }
//
//    if (delta[1] > 0.5 * dd_boxsize[1]) {
//        delta[1] -= dd_boxsize[1];
//    }
//    else if (delta[1] < -0.5 * dd_boxsize[1]) {
//        delta[1] += dd_boxsize[1];
//    }
//
//    if (delta[2] > 0.5 * dd_boxsize[2]) {
//        delta[2] -= dd_boxsize[2];
//    }
//    else if (delta[2] < -0.5 * dd_boxsize[2]) {
//        delta[2] += dd_boxsize[2];
//    }
//
//
//
//    // Calculate the distance between the two particles
//    double r = sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
//
//    // Compute the energy contribution
//    double deltaIdeal = r - params.d;
//    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
//
//    // Accumulate energy to totalPEnergy using atomicAdd
//    atomicAdd(d_totalPEnergy, energy);  // You can switch back to atomicAdd if needed
//
//    // Compute the derivative of the energy with respect to the distance
//    double dEdR = params.k * deltaIdeal;
//
//    // Normalize the delta vector and scale by dEdR
//    Coords3D force;
//    if (r > 0) {
//        force[0] = delta[0] * (dEdR / r);
//        force[1] = delta[1] * (dEdR / r);
//        force[2] = delta[2] * (dEdR / r);
//    }
//    else {
//        force = Coords3D(0.0, 0.0, 0.0);  // Prevent division by zero
//    }
//
//    // Update forces for both atoms using atomicAdd to avoid race conditions
//    //atomicAdd(&forces[params.p1], force);
//    //atomicAdd(&forces[params.p2], -force);
//
//    //atomicAddCoords(forces[params.p1], force);
//    //atomicAddCoords(forces[params.p2], Coords3D{ -force[0], -force[1], -force[2]});
//        // Update forces for both atoms using atomicAdd to avoid race conditions
//
//
//    atomicAdd(&d_forces[params.p1][0], force[0]);
//    atomicAdd(&d_forces[params.p1][1], force[1]);
//    atomicAdd(&d_forces[params.p1][2], force[2]);
//    atomicAdd(&d_forces[params.p2][0], -force[0]);
//    atomicAdd(&d_forces[params.p2][1], -force[1]);
//    atomicAdd(&d_forces[params.p2][2], -force[2]);
//}
//
//
//
//
//void launchKernelBondForces2(Coords3D* atomPositions,HBondParams* bondParams,Coords3D* forces,double* totalPEnergy, Coords3D* boxInfo,int numAtoms, int numBonds) {
//
//    // Device memory allocation
//    Coords3D* d_atomPositions;
//    HBondParams* d_bondParams;
//    Coords3D* d_forces;
//    double* d_totalPEnergy;
//    Coords3D* d_boxsize;
//
//
//    cudaMalloc(&d_atomPositions, numAtoms * sizeof(Coords3D));
//    cudaMalloc(&d_bondParams, numBonds * sizeof(HBondParams));
//    cudaMalloc(&d_forces, numAtoms * sizeof(Coords3D));
//    cudaMalloc(&d_totalPEnergy, sizeof(double));
//    cudaMalloc(&d_boxsize, sizeof(Coords3D));
//
//    // Copy data from host to device
//    cudaMemcpy(d_atomPositions, atomPositions, numAtoms * sizeof(Coords3D), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_bondParams, bondParams, numBonds * sizeof(HBondParams), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_forces, forces, numAtoms * sizeof(Coords3D), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_totalPEnergy, totalPEnergy, sizeof(double), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_boxsize, boxInfo, sizeof(Coords3D), cudaMemcpyHostToDevice);
//
//    // Define block size and number of blocks
//    int blockSize = 32; // Example block size
//    int numBlocks = (numBonds + blockSize - 1) / blockSize;
//
//    // Launch the kernel using float3 types on the device
//    BondForcesKernel_global2 << < numBlocks, blockSize >> > (d_atomPositions, d_bondParams, d_forces, d_totalPEnergy, numBonds, d_boxsize);
//    //BondForcesKernel_global2 << < numBlocks, blockSize >> > (d_atomPositions, d_bondParams,d_forces, d_totalPEnergy, numBonds, d_boxsize);
//
//    // Copy results back to host
//    cudaMemcpy(forces, d_forces, numAtoms * sizeof(Coords3D), cudaMemcpyDeviceToHost);
//    cudaMemcpy(totalPEnergy, d_totalPEnergy, sizeof(double), cudaMemcpyDeviceToHost);
//
//
//    // Free device memory
//    cudaFree(d_atomPositions);
//    cudaFree(d_bondParams);
//    cudaFree(d_forces);
//    cudaFree(d_totalPEnergy);
//    cudaFree(d_boxsize);
//}
