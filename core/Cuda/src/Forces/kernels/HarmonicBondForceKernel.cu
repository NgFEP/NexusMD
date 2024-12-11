// HarmonicBondForce.cu
#include <PeriodicBoundaryConditionKernel.h>
#include "HarmonicBondForceKernel.h"
#include <vector>

int adjustBlockSize(int blockSize, int maxThreadsPerBlock) {
    // Find the closest multiple of 32 above or equal to the current block size
    if (blockSize % 32 != 0) {
        blockSize = ((blockSize / 32) + 1) * 32;
    }

    // Ensure it does not exceed the maximum allowed threads per block
    if (blockSize > maxThreadsPerBlock) {
        blockSize = maxThreadsPerBlock;
    }

    return blockSize;
}

//__global__ void BondForcesKernel_shared(
//    double3* d_atomPositions,
//    BondParams* d_bondParams,
//    double3* d_forces,
//    double* d_totalPEnergy,
//    double3* d_boxsize,
//    ModifiedAtomBondInfo* d_atomsBondLoaded,
//    int numAtomsBondLoaded,
//    int maxBondsPerBlock
//) {
//    // Define the shared memory block
//    extern __shared__ char sharedMemory[];
//
//    // Split the shared memory into separate sections
//    double3* sharedPositions = (double3*)sharedMemory; // Positions of bond atoms
//    BondParams* sharedBondParams = (BondParams*)(sharedPositions + 2 * maxBondsPerBlock);
//    int* sharedAtomIndices = (int*)(sharedBondParams + maxBondsPerBlock);
//
//    // Determine thread and block IDs
//    int tid = threadIdx.x;
//    int bondsLoaded = 0;  // Counter for loaded bonds
//    int atomsLoaded = 0;  // Counter for loaded atoms
//
//    // Initialize shared memory for atom indices
//    if (tid < maxBondsPerBlock * 2) {
//        sharedAtomIndices[tid] = -1; // Initialize to an invalid index
//    }
//    __syncthreads();
//
//    // Loop over atoms until we load enough bonds
//    for (int atomIdx = 0; atomIdx < numAtomsBondLoaded; atomIdx++) {
//        if (bondsLoaded >= maxBondsPerBlock) break;
//
//        // Load atom information
//        ModifiedAtomBondInfo atomInfo = d_atomsBondLoaded[atomIdx];
//
//        // Flag to indicate if the current atom itself is loaded into shared memory
//        bool mainAtomLoaded = false;
//
//        // First, add the main atom if it's not already loaded
//        for (int i = 0; i < atomsLoaded; i++) {
//            if (sharedAtomIndices[i] == atomInfo.atomIdx) {
//                mainAtomLoaded = true;
//                break;
//            }
//        }
//        if (!mainAtomLoaded && tid == atomsLoaded) {
//            sharedPositions[atomsLoaded] = d_atomPositions[atomInfo.atomIdx];
//            sharedAtomIndices[atomsLoaded] = atomInfo.atomIdx;
//            atomsLoaded++;
//        }
//        __syncthreads(); // Ensure the main atom is loaded
//
//        // Load all bonds associated with this atom
//        for (int bondIdx = 0; bondIdx < atomInfo.numBonds && bondsLoaded < maxBondsPerBlock; bondIdx++) {
//            int bondId = atomInfo.bonds[bondIdx];
//
//            // Load positions for the second atom in the bond
//            int atom2 = d_bondParams[bondId].p1 == atomInfo.atomIdx ? d_bondParams[bondId].p2 : d_bondParams[bondId].p1;
//
//            // Check if atom2 is already loaded
//            bool atom2Loaded = false;
//            for (int i = 0; i < atomsLoaded; i++) {
//                if (sharedAtomIndices[i] == atom2) {
//                    atom2Loaded = true;
//                    break;
//                }
//            }
//
//            // If not loaded, load atom2 into shared memory
//            if (!atom2Loaded && tid == atomsLoaded) {
//                sharedPositions[atomsLoaded] = d_atomPositions[atom2];
//                sharedAtomIndices[atomsLoaded] = atom2;
//                atomsLoaded++;
//            }
//            __syncthreads(); // Ensure atom2 is loaded before proceeding
//
//            // Load bond parameters into shared memory
//            if (tid == bondsLoaded) {
//                sharedBondParams[bondsLoaded] = d_bondParams[bondId];
//                sharedPositions[2 * bondsLoaded] = d_atomPositions[atomInfo.atomIdx];
//                sharedPositions[2 * bondsLoaded + 1] = d_atomPositions[atom2];
//            }
//            bondsLoaded++;
//            __syncthreads(); // Ensure bond parameters are fully loaded
//        }
//
//        // Process bonds if max bonds are reached or last atom
//        if (bondsLoaded >= maxBondsPerBlock || atomIdx == numAtomsBondLoaded - 1) {
//            if (tid < bondsLoaded) {
//                double3 pos1 = sharedPositions[2 * tid];
//                double3 pos2 = sharedPositions[2 * tid + 1];
//                BondParams params = sharedBondParams[tid];
//
//                double3 delta;
//                minimumImageVectorDevice(&pos1, &pos2, &delta, d_boxsize);
//                double r = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
//
//                double deltaIdeal = r - params.d;
//                double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
//                atomicAdd(d_totalPEnergy, energy);
//
//                double dEdR = params.k * deltaIdeal;
//                double3 force = make_double3(delta.x * (dEdR / r), delta.y * (dEdR / r), delta.z * (dEdR / r));
//
//                // Apply forces to atoms atom1 and atom2
//                int atom1Idx = sharedAtomIndices[2 * tid];
//                int atom2Idx = sharedAtomIndices[2 * tid + 1];
//                atomicAdd(&d_forces[atom1Idx].x, force.x);
//                atomicAdd(&d_forces[atom1Idx].y, force.y);
//                atomicAdd(&d_forces[atom1Idx].z, force.z);
//                atomicAdd(&d_forces[atom2Idx].x, -force.x);
//                atomicAdd(&d_forces[atom2Idx].y, -force.y);
//                atomicAdd(&d_forces[atom2Idx].z, -force.z);
//            }
//            __syncthreads();
//
//            // Reset counters for next batch
//            bondsLoaded = 0;
//            atomsLoaded = 0;
//
//            // Reinitialize atom indices for the next batch
//            if (tid < maxBondsPerBlock * 2) {
//                sharedAtomIndices[tid] = -1;
//            }
//            __syncthreads();
//        }
//    }
//
//
//
//    //// Determine thread and block IDs
//    //int tid = threadIdx.x;
//    //int bondsLoaded = 0; // Counter for loaded bonds
//
//    //// Loop over atoms assigned to this block
//    //for (int atomIdx = blockIdx.x; atomIdx < numAtomsBondLoaded; atomIdx += gridDim.x) {
//    //    ModifiedAtomBondInfo atomInfo = d_atomsBondLoaded[atomIdx];
//
//    //    // Load atom information into shared memory
//    //    if (tid == 0) {
//    //        sharedAtomsInfo[0] = atomInfo;
//    //    }
//    //    __syncthreads();
//
//    //    // Load all bonds associated with this atom based on numBonds
//    //    for (int bondIdx = 0; bondIdx < atomInfo.numBonds && bondsLoaded < maxBondsPerBlock; bondIdx++) {
//    //        int bondId = atomInfo.bonds[bondIdx];
//
//    //        // Ensure the bond index is valid (within bounds)
//    //        if (bondIdx < atomInfo.numBonds) {
//    //            if (tid == bondsLoaded) {
//    //                sharedPositions[2 * bondsLoaded] = d_atomPositions[d_bondParams[bondId].p1];
//    //                sharedPositions[2 * bondsLoaded + 1] = d_atomPositions[d_bondParams[bondId].p2];
//    //                sharedBondParams[bondsLoaded] = d_bondParams[bondId];
//    //            }
//    //            bondsLoaded++;
//    //        }
//    //    }
//    //    __syncthreads();
//
//    //    // Each thread processes one bond
//    //    if (tid < bondsLoaded) {
//    //        double3 pos1 = sharedPositions[2 * tid];
//    //        double3 pos2 = sharedPositions[2 * tid + 1];
//    //        BondParams params = sharedBondParams[tid];
//
//    //        double3 delta;
//    //        minimumImageVectorDevice(&pos1, &pos2, &delta, d_boxsize);
//    //        double r = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
//
//    //        double deltaIdeal = r - params.d;
//    //        double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
//    //        atomicAdd(d_totalPEnergy, energy);
//
//    //        double dEdR = params.k * deltaIdeal;
//    //        double3 force = make_double3(delta.x * (dEdR / r), delta.y * (dEdR / r), delta.z * (dEdR / r));
//
//    //        // Apply forces to atoms
//    //        atomicAdd(&d_forces[d_bondParams[tid].p1].x, force.x);
//    //        atomicAdd(&d_forces[d_bondParams[tid].p1].y, force.y);
//    //        atomicAdd(&d_forces[d_bondParams[tid].p1].z, force.z);
//    //        atomicAdd(&d_forces[d_bondParams[tid].p2].x, -force.x);
//    //        atomicAdd(&d_forces[d_bondParams[tid].p2].y, -force.y);
//    //        atomicAdd(&d_forces[d_bondParams[tid].p2].z, -force.z);
//    //    }
//    //    __syncthreads(); // Ensure all threads finish before starting the next atom
//
//    //    bondsLoaded = 0; // Reset bond counter for the next atom
//    //}
//}
//
//
//
//
//void launchKernelBondForcesShared(
//    double3* d_atomPositions,
//    BondParams* d_bondParams,
//    double3* d_forces,
//    double* d_totalPEnergy,
//    double3* d_boxsize,
//    D_Residues* d_desidues,
//    int numResidues
//) {
//    // Calculate max bonds per block based on shared memory availability
//    const int sizeOfPositions = sizeof(double3);
//    const int sizeOfBondParams = sizeof(BondParams);
//    const int sizeOfResidue = sizeof(D_Residues);
//
//    int deviceId;
//    cudaGetDevice(&deviceId); // Get the current device ID
//
//    // Automatically obtain shared memory per block from the GPU
//    int sharedMemPerBlock;
//    cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);
//
//    //int deviceId = 0;
//    //cudaError_t status;
//
//    //// Get the current device ID
//    //status = cudaGetDevice(&deviceId);
//    //if (status != cudaSuccess) {
//    //    fprintf(stderr, "Failed to get current device ID: %s\n", cudaGetErrorString(status));
//    //    deviceId = 0; // Default to device 0 if there's an error
//    //}
//
//    //// Automatically obtain shared memory per block from the GPU
//    //int sharedMemPerBlock = 0;
//    //status = cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);
//    //if (status != cudaSuccess) {
//    //    fprintf(stderr, "Failed to get max shared memory per block: %s\n", cudaGetErrorString(status));
//    //    sharedMemPerBlock = 49152; // Default to 48 KB if the query fails
//    //}
//
//    // Print out the value for debugging purposes
//    printf("Shared memory per block for device %d: %d bytes\n", deviceId, sharedMemPerBlock);
//
//
//    // Calculate maxBondsPerBlock with a 10% margin added
//    int maxBondsPerBlock = (int)(sharedMemPerBlock / (2 * sizeOfPositions + sizeOfBondParams + sizeOfResidue) * 1.1); // Adding 10% margin
//    int blockSize = maxBondsPerBlock; // Adaptive block size based on memory availability
//
//
//
//    // Calculate shared memory size needed with a 10% margin
//    int sharedMemorySize = (int)((2 * sizeOfPositions * maxBondsPerBlock) +
//        (sizeOfBondParams * maxBondsPerBlock) +
//        sizeOfResidue * 1.1); // Adding 10% margin to the shared memory size
//
//    // Calculate number of blocks
//    int numBlocks = (numResidues + blockSize - 1) / blockSize;
//
//
//    int maxThreadsPerBlock;
//    cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, deviceId);
//
//    // Adjusting blockSize
//    blockSize = adjustBlockSize(blockSize, maxThreadsPerBlock);
//    // Recalculating parameters
//    maxBondsPerBlock = blockSize;
//    // Calculate shared memory size needed with a 10% margin
//    sharedMemorySize = (int)((2 * sizeOfPositions * maxBondsPerBlock) +
//        (sizeOfBondParams * maxBondsPerBlock) +
//        sizeOfResidue * 1.1);
//
//    numBlocks = (numResidues + blockSize - 1) / blockSize;
//
//
//    printf("maxBondsPerBlock %d, blockSize %d,maxThreadsPerBlock %d, sharedMemorySize %d, numBlocks %d\n", maxBondsPerBlock, blockSize, maxThreadsPerBlock, sharedMemorySize, numBlocks);
//
//    // Launch the kernel
//    BondForcesKernel_shared << <numBlocks, blockSize, sharedMemorySize >> > (
//        d_atomPositions,
//        d_bondParams,
//        d_forces,
//        d_totalPEnergy,
//        d_boxsize,
//        d_residues,
//        numResidues,
//        maxBondsPerBlock
//        );
//
//    // Check for errors
//    cudaError_t err = cudaGetLastError();
//    if (err != cudaSuccess) {
//        printf("Error launching BondForcesKernel_shared: %s\n", cudaGetErrorString(err));
//    }
//}
//int calculateResidueMemory(const Residues& residue) {
//    // Calculate dynamic array sizes
//    int atomMemory = residue.AllAtomsCount * sizeof(int);
//    int bondMemory = residue.AllBondsCount * sizeof(int);
//    int hAtomMemory = residue.HAtomsCount * sizeof(int);
//    int nonHAtomMemory = residue.NonHAtomsCount * sizeof(int);
//    int hBondMemory = residue.HBondsCount * sizeof(int);
//    int nonHBondMemory = residue.NonHBondsCount * sizeof(int);
//    int nonResBondAtomMemory = residue.NonResBondAtomsCount * sizeof(int);
//
//    // Pointer and static member contributions
//    int staticMemory = sizeof(residue.resName);
//
//    // Total memory for this residue
//    int totalMemory = atomMemory + bondMemory + hAtomMemory +
//        nonHAtomMemory + hBondMemory + nonHBondMemory +
//        nonResBondAtomMemory + staticMemory;
//
//    return totalMemory;
//}


//__global__ void BondForcesKernel_shared(double3* d_atomPositions,
//    BondParams* d_bondParams,
//    double3* d_forces,
//    double* d_totalPEnergy,
//    double3* d_boxsize,
//    D_Residues* d_residues, 
//    int sharedMemorySize,
//    int* residueBondCounts, 
//    int numResidues, 
//    int maxBondsPerBlock) 
//{
//    extern __shared__ D_Residues sharedResidues[]; // Shared memory for residues
//
//    int blockId = blockIdx.x;
//    int threadId = threadIdx.x;
//
//    // Step 1: Identify range of residues for this block
//    int startResidue = 0;
//    int occupiedMemory = 0;
//
//    // Locate the starting residue for this block
//    for (int i = 0; i < numResidues; ++i) {
//        if (occupiedMemory >= sharedMemorySize) {//blockId * maxBondsPerBlock
//            startResidue = i;
//            break;
//        }
//        occupiedMemory += d_residues[i].resMemSize;
//    }
//
//    // Accumulate residues for this block while respecting maxBondsPerBlock
//    int endResidue = startResidue;
//    occupiedMemory = 0; // Reset for the block
//    while (endResidue < numResidues && occupiedMemory + residueBondCounts[endResidue] <= maxBondsPerBlock) {
//        occupiedMemory += residueBondCounts[endResidue];
//        ++endResidue;
//    }
//
//    // Total number of residues for this block
//    int blockResidueCount = endResidue - startResidue;
//
//    // Step 2: Load residues into shared memory (coalesced access)
//    for (int i = threadId; i < blockResidueCount; i += blockDim.x) {
//        sharedResidues[i] = d_residues[startResidue + i];
//    }
//
//    __syncthreads(); // Ensure all data is loaded
//
//    // Step 3: Process residues in shared memory
//    if (threadId < blockResidueCount) {
//        D_Residues r = sharedResidues[threadId];
//
//        if (r.type == 0) {
//            // Protein residue processing
//            for (int b = 0; b < r.bondCount; b++) {
//                int bondIndex = r.bondIndices[b];
//                // Perform protein-specific computations...
//            }
//        }
//        else if (r.type == 1) {
//            // Water molecule processing
//            for (int b = 0; b < r.bondCount; b++) {
//                int bondIndex = r.bondIndices[b];
//                // Perform water-specific computations...
//            }
//        }
//    }
//
//    __syncthreads(); // Optional: Synchronize before the next block-level operation
//}

struct AtomPosition {
    int globalAtomIdx;  // Global atom index
    double3 position;   // Position of the atom
};


__global__ void BondForcesKernel_shared(
    double3* d_atomPositions,
    BondParams* d_bondParams,
    double3* d_forces,
    double* d_totalPEnergy,
    double3* d_boxsize,
    D_Residues* d_residues,
    int* d_startResidues,
    int* d_endResidues
) {
    extern __shared__ char sharedMemory[];  // Dynamically allocated shared memory
    //printf("atom1Pos");
    int blockId = blockIdx.x;
    int threadId = threadIdx.x;

    int startResidue = d_startResidues[blockId];
    int endResidue = d_endResidues[blockId];

    // Shared memory pointers
    AtomPosition* sharedAtoms = (AtomPosition*)sharedMemory;


    for (int resIdx = startResidue; resIdx <= endResidue; ++resIdx) {
        D_Residues residue = d_residues[resIdx];

        // Copy atom positions and indices for this residue into shared memory
        for (int i = threadId; i < residue.AllAtomsCount; i += blockDim.x) {
            int globalAtomIdx = residue.AllAtomsIndices[i];
            sharedAtoms[i].globalAtomIdx = globalAtomIdx;
            sharedAtoms[i].position = d_atomPositions[globalAtomIdx];
        }

        __syncthreads();  // Ensure all threads have loaded data

        // Calculate forces for each bond
        for (int bondIdx = threadId; bondIdx < residue.AllBondsCount; bondIdx += blockDim.x) {
            int bondGlobalIdx = residue.AllBondsIndices[bondIdx];
            BondParams bond = d_bondParams[bondGlobalIdx];

            int atom1Idx = bond.p1;
            int atom2Idx = bond.p2;

            // Simulating dictionary access by searching for atom1Idx and atom2Idx in shared memory
            double3 atom1Pos; // = { 0.0, 0.0, 0.0 };
            double3 atom2Pos; // = { 0.0, 0.0, 0.0 };

            // Search for atom1 position in shared memory
            for (int i = 0; i < residue.AllAtomsCount; ++i) {
                if (sharedAtoms[i].globalAtomIdx == atom1Idx) {
                    atom1Pos = sharedAtoms[i].position;
                    break;
                }
            }

            //printf("atom1Pos.x: %f, atom1Pos.y: %f, atom1Pos.z: %f\n", atom1Pos.x, atom1Pos.y, atom1Pos.y);


            // Search for atom2 position in shared memory
            for (int i = 0; i < residue.AllAtomsCount; ++i) {
                if (sharedAtoms[i].globalAtomIdx == atom2Idx) {
                    atom2Pos = sharedAtoms[i].position;
                    break;
                }
            }

            //printf("atom2Pos.x: %f, atom2Pos.y: %f, atom2Pos.z: %f\n", atom2Pos.x, atom2Pos.y, atom2Pos.y);


            // Compute forces based on atom positions
            double3 delta;
            minimumImageVectorDevice(&atom1Pos, &atom2Pos, &delta, d_boxsize);

            // Calculate the distance between the two particles
            double r = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);

            // Compute the energy contribution
            double deltaIdeal = r - bond.d;
            double energy = 0.5 * bond.k * deltaIdeal * deltaIdeal;



            // Compute the derivative of the energy with respect to the distance
            double dEdR = bond.k * deltaIdeal;

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

            if (residue.resName == "WAT") {
                atomicAdd(&d_forces[atom1Idx].x, force.x);
                atomicAdd(&d_forces[atom1Idx].y, force.y);
                atomicAdd(&d_forces[atom1Idx].z, force.z);
                atomicAdd(&d_forces[atom2Idx].x, -force.x);
                atomicAdd(&d_forces[atom2Idx].y, -force.y);
                atomicAdd(&d_forces[atom2Idx].z, -force.z);
                // Accumulate energy to totalPEnergy using atomicAdd
                atomicAdd(d_totalPEnergy, energy);  // You can switch back to atomicAdd if needed
            }
            else {
                int counter = 0;
                for (int i = 0; i < residue.NonResBondAtomsCount; ++i) {
                    if (residue.NonResBondAtoms[i] == atom2Idx) {
                        atomicAdd(&d_forces[atom1Idx].x, force.x);
                        atomicAdd(&d_forces[atom1Idx].y, force.y);
                        atomicAdd(&d_forces[atom1Idx].z, force.z);
                        atomicAdd(d_totalPEnergy, 0.5 * energy); 
                        break;
                    }
                    else if (residue.NonResBondAtoms[i] == atom1Idx) {
                        atomicAdd(&d_forces[atom2Idx].x, -force.x);
                        atomicAdd(&d_forces[atom2Idx].y, -force.y);
                        atomicAdd(&d_forces[atom2Idx].z, -force.z);
                        atomicAdd(d_totalPEnergy, 0.5 * energy);  
                        break;
                    }
                    counter++;
                }
                if (counter == residue.NonResBondAtomsCount)
                {
                    atomicAdd(&d_forces[atom1Idx].x, force.x);
                    atomicAdd(&d_forces[atom1Idx].y, force.y);
                    atomicAdd(&d_forces[atom1Idx].z, force.z);
                    atomicAdd(&d_forces[atom2Idx].x, -force.x);
                    atomicAdd(&d_forces[atom2Idx].y, -force.y);
                    atomicAdd(&d_forces[atom2Idx].z, -force.z);
                    // Accumulate energy to totalPEnergy using atomicAdd
                    atomicAdd(d_totalPEnergy, energy);  // You can switch back to atomicAdd if needed
                }
            }


            //// Accumulate potential energy
            //if (threadId == 0) {
            //    atomicAdd(d_totalPEnergy, 0.5 * bond.k * deltaIdeal * deltaIdeal);
            //}
        }

        __syncthreads();
    }
}



void launchKernelBondForcesShared(
    double3* d_atomPositions,
    BondParams* d_bondParams,
    double3* d_forces,
    double* d_totalPEnergy,
    double3* d_boxsize,
    D_Residues* d_residues,
    int* d_startResidues,          // Start indices for residues in each block
    int* d_endResidues,            // End indices for residues in each block
    int _numBlocks,
    int totalBondsInResidues

) {
    // Calculate sizes of static components
    const int sizeOfPositions = sizeof(double3);
    const int sizeOfBondParams = sizeof(BondParams);
    //const int maxResidueSize= sizeof(h_residue);
    // Calculate the maximum memory required for a single residue
    //int maxResidueSize = calculateResidueMemory(h_residue); // Use the first residue (largest)

    // Get current device ID
    int deviceId;
    cudaGetDevice(&deviceId);

    //// Obtain shared memory per block from the GPU
    int sharedMemPerBlock;
    cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);

    // Obtain max threads per block for the device
    int maxThreadsPerBlock;
    cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, deviceId);

    //// Calculate maximum memory per bond
    //int maxMemoryPerBond = 2 * sizeOfPositions + sizeOfBondParams + maxResidueSize;

    //// Calculate max bonds per block with a margin (no more than maxThreadsPerBlock)
    //int maxBondsPerBlock = sharedMemPerBlock / maxMemoryPerBond;
    //maxBondsPerBlock = std::min(maxBondsPerBlock, maxThreadsPerBlock);

    //// Adaptive block size based on available resources
    //int blockSize = maxBondsPerBlock;

    //// Calculate shared memory size for the kernel
    //int sharedMemorySize = blockSize * maxMemoryPerBond;
    //// Calculate number of blocks required
    //int numBlocks = (totalBondsInResidues + blockSize - 1) / blockSize;


    //int totalMemoryAllBonds = totalBondsInResidues * (2 * sizeOfPositions + sizeOfBondParams) + totalResiduesSize;
    //int numBlocks = 1.1* totalMemoryAllBonds / sharedMemPerBlock; // 10% margin

    ////numBlocks = (totalBondsInResidues + blockSize - 1) / blockSize;
    int sharedMemorySize = 0.95* sharedMemPerBlock;
    int numBlocks = _numBlocks;

    int blockSize = totalBondsInResidues / numBlocks;

    blockSize = adjustBlockSize(blockSize, maxThreadsPerBlock);





    // Debugging information
    //printf("Device ID: %d\n", deviceId);
    //printf("Shared memory per block: %d bytes\n", sharedMemPerBlock);
    //printf("Max threads per block: %d\n", maxThreadsPerBlock);
    //printf("Max bonds per block: %d\n", maxBondsPerBlock);
    //printf("Block size: %d\n", blockSize);
    //printf("Shared memory size: %d bytes\n", sharedMemorySize);
    //printf("Number of blocks: %d\n", numBlocks);

    // Launch the kernel
    BondForcesKernel_shared << <numBlocks, blockSize, sharedMemorySize >> > (
        d_atomPositions,
        d_bondParams,
        d_forces,
        d_totalPEnergy,
        d_boxsize,
        d_residues,
        d_startResidues,
        d_endResidues
        );

    // Check for errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Error launching BondForcesKernel_shared: %s\n", cudaGetErrorString(err));
    }
}


//void launchKernelBondForcesShared(
//    double3* d_atomPositions,
//    BondParams* d_bondParams,
//    double3* d_forces,
//    double* d_totalPEnergy,
//    double3* d_boxsize,
//    PResidues* d_pResidues,
//    WResidues* d_wResidues,
//    int numPResidues,
//    int numWResidues
//) {
//    // Sizes of main components
//    const int sizeOfPositions = sizeof(double3);
//    const int sizeOfBondParams = sizeof(BondParams);
//    const int sizeOfPResidue = sizeof(PResidues);
//    const int sizeOfWResidue = sizeof(WResidues);
//
//    // Get device properties
//    int deviceId;
//    cudaGetDevice(&deviceId);
//
//    int sharedMemPerBlock;
//    cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);
//
//    printf("Shared memory per block for device %d: %d bytes\n", deviceId, sharedMemPerBlock);
//
//    // Calculate memory requirements for PResidues and WResidues
//    int pResidueMem = numPResidues * sizeOfPResidue;
//    int wResidueMem = numWResidues * sizeOfWResidue;
//
//    // Add contributions from optional fields (if required)
//    // Adjust this if optional fields are stored in global memory or accessed differently
//    int maxOptionalMem = 0;
//    maxOptionalMem += sizeof(int) * numPResidues * 6; // Assuming optional vectors with integers in PResidues
//    maxOptionalMem += sizeof(int) * numWResidues * 3; // Assuming optional vectors with integers in WResidues
//
//    // Total shared memory needed per block
//    int sharedMemorySize = (2 * sizeOfPositions) + sizeOfBondParams + pResidueMem + wResidueMem + maxOptionalMem;
//
//    if (sharedMemorySize > sharedMemPerBlock) {
//        printf("Error: Shared memory requirement (%d bytes) exceeds available memory (%d bytes)\n",
//            sharedMemorySize, sharedMemPerBlock);
//        return;
//    }
//
//    // Calculate max bonds per block
//    int maxBondsPerBlock = sharedMemPerBlock / sharedMemorySize;
//
//    // Determine block and grid size
//    int maxThreadsPerBlock;
//    cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, deviceId);
//
//    int blockSize = adjustBlockSize(maxBondsPerBlock, maxThreadsPerBlock);
//    int numBlocks = (numPResidues + numWResidues + blockSize - 1) / blockSize;
//
//    printf("maxBondsPerBlock: %d, blockSize: %d, numBlocks: %d, sharedMemorySize: %d\n",
//        maxBondsPerBlock, blockSize, numBlocks, sharedMemorySize);
//
//    // Launch the kernel
//    //BondForcesKernel_shared << <numBlocks, blockSize, sharedMemorySize >> > (
//    //    d_atomPositions,
//    //    d_bondParams,
//    //    d_forces,
//    //    d_totalPEnergy,
//    //    d_boxsize,
//    //    d_pResidues,
//    //    d_wResidues,
//    //    numPResidues,
//    //    numWResidues,
//    //    maxBondsPerBlock
//    //    );
//
//    // Check for errors
//    cudaError_t err = cudaGetLastError();
//    if (err != cudaSuccess) {
//        printf("Error launching BondForcesKernel_shared: %s\n", cudaGetErrorString(err));
//    }
//}









__global__ void BondForcesKernel_global(
    double3* d_atomPositions,
    BondParams* d_bondParams,
    double3* d_forces,
    double* d_totalPEnergy,
    int d_numBonds,
    double3* d_boxsize)  // Periodic boundary conditions info
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= d_numBonds) return;  // Ensure we don't access out-of-bounds

    // Load bond parameters from global memory
    BondParams& params = d_bondParams[idx];

    // Load atom positions from global memory
    double3 coords1 = d_atomPositions[params.p1];
    double3 coords2 = d_atomPositions[params.p2];
    //printf("params.p1: %f, params.p2: %f\n", params.p2, params.p2);
    //printf("coords1.x: %f, coords1.y: %f, coords1.z: %f\n", coords1.x, coords1.y, coords1.z);

    // Compute the vector from pos1 to pos2 considering periodic boundary conditions
    double3 delta;
    minimumImageVectorDevice(&coords1, &coords2, &delta, d_boxsize);

    // Calculate the distance between the two particles
    double r = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);

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
void launchKernelBondForcesGlobal(
    double3* d_atomPositions,  // Host-side Coords3D
    BondParams* d_bondParams,
    double3* d_forces,
    double* d_totalPEnergy,
    double3* d_boxsize,
    int _numBonds
) {

    // Define block size and number of blocks
    int blockSize = 256; // Example block size
    int numBlocks = (_numBonds + blockSize - 1) / blockSize;

    // Launch the kernel using double3 types on the device
    BondForcesKernel_global << < numBlocks, blockSize >> > (d_atomPositions, d_bondParams, d_forces, d_totalPEnergy, _numBonds, d_boxsize);

    //cudaError_t error = cudaGetLastError();
    //if (error != cudaSuccess) {
    //    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
    //}
    //cudaDeviceSynchronize(); // Wait for the kernel to finish

}








//. using a  __device__ const int blocksize = 320; Is faster than using blockDim.x all the time
//__device__ const int blockSize = 256;

//__global__ void BondForcesKernel_shared(
//    const Coords3D* atomPositions,
//    const BondParams* bondParams,
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
//    __shared__ BondParams sharedBondParams[blockSize];  // Assuming max 256 threads/block
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
//    const BondParams& params = sharedBondParams[threadIdx.x];
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



//__global__ void BondForcesKernel_global(
//    double3* d_atomPositions,
//    BondParams* d_bondParams,
//    double3* d_forces,
//    double* d_totalPEnergy,
//    int d_numBonds,
//    double3* d_boxsize)  // Periodic boundary conditions info
//{
//    int idx = threadIdx.x + blockIdx.x * blockDim.x;
//
//    if (idx >= d_numBonds) return;  // Ensure we don't access out-of-bounds
//
//    // Load bond parameters from global memory
//    BondParams& params = d_bondParams[idx];
//
//    // Load atom positions from global memory
//    double3 coords1 = d_atomPositions[params.p1];
//    double3 coords2 = d_atomPositions[params.p2];
//    //printf("params.p1: %f, params.p2: %f\n", params.p2, params.p2);
//    //printf("coords1.x: %f, coords1.y: %f, coords1.z: %f\n", coords1.x, coords1.y, coords1.z);
//
//    // Compute the vector from pos1 to pos2 considering periodic boundary conditions
//    double3 delta;
//    minimumImageVectorDevice(&coords1, &coords2, &delta, d_boxsize);
//
//    // Calculate the distance between the two particles
//    double r = sqrt(delta.x*delta.x + delta.y * delta.y + delta.z * delta.z);
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
//    double3 force;
//    if (r > 0) {
//        force.x = delta.x * (dEdR / r);
//        force.y = delta.y * (dEdR / r);
//        force.z = delta.z * (dEdR / r);
//    }
//    else {
//        force = make_double3(0.0, 0.0, 0.0);  // Prevent division by zero
//    }
//
//    // Update forces for both atoms using atomicAdd to avoid race conditions
//    //atomicAdd(&forces[params.p1], force);
//    //atomicAdd(&forces[params.p2], -force);
//
//    //atomicAddCoords(forces[params.p1], force);
//    //atomicAddCoords(forces[params.p2], Coords3D{ -force[0], -force[1], -force[2]});
//    // Update forces for both atoms using atomicAdd to avoid race conditions
//
//
//    atomicAdd(&d_forces[params.p1].x, force.x);
//    atomicAdd(&d_forces[params.p1].y, force.y);
//    atomicAdd(&d_forces[params.p1].z, force.z);
//    atomicAdd(&d_forces[params.p2].x, -force.x);
//    atomicAdd(&d_forces[params.p2].y, -force.y);
//    atomicAdd(&d_forces[params.p2].z, -force.z);
//}
//
//
////correct for double
//void launchKernelBondForces(
//    double3* atomPositions,  // Host-side Coords3D
//    BondParams* bondParams,
//    double3* forces,
//    double* totalPEnergy,
//    double3* boxInfo,
//    int numAtoms,
//    int numBonds
//) {
//
//    // Device memory allocation
//    double3* d_atomPositions;
//    BondParams* d_bondParams;
//    double3* d_forces;
//    double* d_totalPEnergy;
//    double3* d_boxsize;
//
//
//    cudaMalloc(&d_atomPositions, numAtoms * sizeof(double3));
//    cudaMalloc(&d_bondParams, numBonds * sizeof(BondParams));
//    cudaMalloc(&d_forces, numAtoms * sizeof(double3));
//    cudaMalloc(&d_totalPEnergy, sizeof(double));
//    cudaMalloc(&d_boxsize, sizeof(double3));
//
//    // Copy data from host to device
//    cudaMemcpy(d_atomPositions, atomPositions, numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_bondParams, bondParams, numBonds * sizeof(BondParams), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_forces, forces, numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_totalPEnergy, totalPEnergy, sizeof(double), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_boxsize, boxInfo, sizeof(double3), cudaMemcpyHostToDevice);
//
//    // Define block size and number of blocks
//    int blockSize = 256; // Example block size
//    int numBlocks = (numBonds + blockSize - 1) / blockSize;
//
//    // Launch the kernel using double3 types on the device
//    BondForcesKernel_global <<< numBlocks, blockSize >>> (d_atomPositions, d_bondParams, d_forces, d_totalPEnergy, numBonds, d_boxsize);
//    //BondForcesKernel_global2 << < numBlocks, blockSize >> > (d_atomPositions, d_bondParams,d_forces, d_totalPEnergy, numBonds, d_boxsize);
//
//    // Copy results back to host
//    cudaMemcpy(forces, d_forces, numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
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
//











//__global__ void BondForcesKernel_global2(
//    Coords3D* d_atomPositions,
//    BondParams* d_bondParams,
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
//    BondParams& params = d_bondParams[idx];
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
//void launchKernelBondForces2(Coords3D* atomPositions,BondParams* bondParams,Coords3D* forces,double* totalPEnergy, Coords3D* boxInfo,int numAtoms, int numBonds) {
//
//    // Device memory allocation
//    Coords3D* d_atomPositions;
//    BondParams* d_bondParams;
//    Coords3D* d_forces;
//    double* d_totalPEnergy;
//    Coords3D* d_boxsize;
//
//
//    cudaMalloc(&d_atomPositions, numAtoms * sizeof(Coords3D));
//    cudaMalloc(&d_bondParams, numBonds * sizeof(BondParams));
//    cudaMalloc(&d_forces, numAtoms * sizeof(Coords3D));
//    cudaMalloc(&d_totalPEnergy, sizeof(double));
//    cudaMalloc(&d_boxsize, sizeof(Coords3D));
//
//    // Copy data from host to device
//    cudaMemcpy(d_atomPositions, atomPositions, numAtoms * sizeof(Coords3D), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_bondParams, bondParams, numBonds * sizeof(BondParams), cudaMemcpyHostToDevice);
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



