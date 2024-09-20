// HarmonicBondForce.cu
#include "Coords3D.h"

//#include <math.h>
#include <stdio.h> // For printf
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <iomanip>//to print with 16 decimal  
#include "PeriodicBoundaryCondition.h"

#include "SystemXMLParser.h"
#include <cuda.h>
#include <cuda_runtime.h> // Stops underlining of __global__
#include <device_launch_parameters.h> // Stops underlining of threadIdx etc.
#include <cooperative_groups.h>


void launchKernelBondForces(double3* atomPositions,  HBondParams* bondParams,double3* forces,double* totalPEnergy,double3* boxInfo, int numAtoms, int numBonds);
//void launchKernelBondForces2(Coords3D* atomPositions, HBondParams* bondParams, Coords3D* forces, double* totalPEnergy, Coords3D* boxInfo, int numAtoms, int numBonds);


//__device__ void minimumImageVector(const Coords3D* pos1, const Coords3D* pos2, Coords3D* delta, const PeriodicBoundaryCondition::BoxInfo* boxInfo);
