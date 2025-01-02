
//#include <math.h>
#include <stdio.h> // For printf
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <algorithm> // For std::min
#include <iomanip>//to print with 16 decimal  
#include "PeriodicBoundaryCondition.h"
#include "SystemXMLParser.h"
#include <cuda.h>
#include <cuda_runtime.h> // Stops underlining of __global__
#include <device_launch_parameters.h> // Stops underlining of threadIdx etc.
#include <cooperative_groups.h>
#include "CudaDataStructures.h"
//#include "PDBResidueParser.h" // Include for residues
//struct ResBondInfo;


//void launchKernelBondForces(double3* atomPositions,  BondParams* bondParams,double3* forces,double* totalPEnergy,double3* boxInfo, int numAtoms, int numBonds);
//void launchKernelBondForces2(Coords3D* atomPositions, BondParams* bondParams, Coords3D* forces, double* totalPEnergy, Coords3D* boxInfo, int numAtoms, int numBonds);


//__device__ void minimumImageVector(const Coords3D* pos1, const Coords3D* pos2, Coords3D* delta, const PeriodicBoundaryCondition::BoxInfo* boxInfo);
void launchKernelBondForcesGlobal(double3* d_atomPositions, BondParams* d_bondParams, double3* d_forces, double* d_totalPEnergy, double3* d_boxsize, int _numBonds);
//void launchKernelBondForcesShared(double3* d_atomPositions, BondParams* d_bondParams, double3* d_forces, double* d_bondPEnergies, double* d_totalPEnergy, double3* d_boxsize, D_Residues* d_residues, int* d_startResidues, int* d_endResidues, int _numBlocks, int totalBondsInResidues);
void launchKernelBondForcesShared(double3* d_atomPositions, double3* d_forces, double* d_totalPEnergy, double3* d_boxsize, D_CudaBonds* d_cudaBonds, int _numBlocks, int totalBondsInResidues);

//void launchKernelBondForcesShared(double3* d_atomPositions, BondParams* d_bondParams, double3* d_forces, double* d_totalPEnergy, double3* d_boxsize, ModifiedAtomBondInfo* d_atomsBondLoaded, int numAtomsBondLoaded);