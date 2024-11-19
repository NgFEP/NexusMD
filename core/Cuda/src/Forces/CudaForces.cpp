//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "CudaForces.h"
#include "utility"

using namespace std;
using namespace Cuda;


void Forces::AddPTorsion(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<PTorsionParams>& torsionParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    //vector<Coords3D> totalForces(atomPositions.size(), Coords3D(0, 0, 0)); // Initialize totalForces with the size of atomPositions

    for (const auto& torsion : torsionParams) {
        //AP is a vector of 4 atoms positions involved in the torsion 
        vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
                                        Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
                                        Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
                                        Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };
        auto forces = PeriodicTorsionForce::calculateForces(AP, torsion, totalPEnergy, boxInfo);
        totalForces[torsion.p1] += forces[0];
        totalForces[torsion.p2] += forces[1];
        totalForces[torsion.p3] += forces[2];
        totalForces[torsion.p4] += forces[3];

        //if (abs(totalForces[torsion.p1][1]) > 200) {
        //    cout << "";
        //}
        //if (abs(totalForces[torsion.p2][1]) > 200) {
        //    cout << "";
        //}
        //if (abs(totalForces[torsion.p3][1]) > 200) {
        //    cout << "";
        //}
        //if (abs(totalForces[torsion.p4][1]) > 200) {
        //    cout << "";
        //}





    }



    ////For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0]
    ////reason for error for (const auto& torsion : torsionParams[0]) { 
    ////torsionParams[0] refers to a single instance of TorsionParameters, not a collection that I can iterate over with a range - based for loop.
    //if (!torsionParams.empty()) {
    //    const auto& torsion = torsionParams[0];

    //    //AP is a vector of 4 atoms positions involved in the torsion 
    //    vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
    //                                 Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
    //                                 Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
    //                                 Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };

    //    auto forces = PeriodicTorsionForce::calculateForces(AP, torsion);
    //    // Assuming forcesAndAngle.second is the vector<Coords3D> of forces
    //    totalForces[torsion.p1] += forces[0];
    //    totalForces[torsion.p2] += forces[1];
    //    totalForces[torsion.p3] += forces[2];
    //    totalForces[torsion.p4] += forces[3];
    //}

}



void Forces::AddHBond(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<BondParams>& bondParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    


    //// No need for manual conversion from Coords3D to double3, pass data directly
    //Coords3D* atomPositions_ptr = const_cast<Coords3D*>(atomPositions.data());
    //Coords3D* forces_ptr = totalForces.data();
    //BondParams* bondParams_ptr = const_cast<BondParams*>(bondParams.data());

    //// Convert boxInfo into Coords3D for the kernel
    //Coords3D boxSize = { boxInfo.boxSize[0], boxInfo.boxSize[1], boxInfo.boxSize[2] };


    //launchKernelBondForces2(atomPositions_ptr, bondParams_ptr, forces_ptr, &totalPEnergy, &boxSize, atomPositions.size(), bondParams.size());


    // Convert Coords3D (host) to double3 (device)
    double3* atomPositions_double3 = new double3[atomPositions.size()];
    double3* forces_double3 = new double3[totalForces.size()];
    BondParams* bondParams_ptr = new BondParams[bondParams.size()];
    double* totalPEnergy_double = new double;
    double3* h_boxSize_double3 = new double3;
    *h_boxSize_double3 = make_double3(boxInfo.boxSize[0], boxInfo.boxSize[1], boxInfo.boxSize[2]);
    *totalPEnergy_double = double(totalPEnergy);


    for (int i = 0; i < atomPositions.size(); ++i) {
        atomPositions_double3[i] = make_double3(atomPositions[i][0], atomPositions[i][1], atomPositions[i][2]);
        forces_double3[i] = make_double3(totalForces[i][0], totalForces[i][1], totalForces[i][2]);
    }
    for (int j = 0; j < bondParams.size(); ++j) {
        bondParams_ptr[j] = bondParams[j];
    }

    //launchKernelBondForces(atomPositions_double3, bondParams_ptr, forces_double3, totalPEnergy_double, h_boxSize_double3, atomPositions.size(), bondParams.size());

    // Convert double3 results back to Coords3D
    for (size_t i = 0; i < atomPositions.size(); ++i) {
        totalForces[i] = Coords3D{ forces_double3[i].x, forces_double3[i].y, forces_double3[i].z };
    }

    totalPEnergy = double(*totalPEnergy_double);


    // Free allocated memory
    delete[] atomPositions_double3;
    delete[] forces_double3;
    delete[] bondParams_ptr;
    delete totalPEnergy_double;
    delete h_boxSize_double3;

    //// Use smart pointers for memory management
    //std::unique_ptr<double3[]> atomPositions_double3(new double3[atomPositions.size()]);
    //std::unique_ptr<double3[]> forces_double3(new double3[totalForces.size()]);
    //std::unique_ptr<
    // []> bondParams_ptr(new BondParams[bondParams.size()]);
    //std::unique_ptr<double> totalPEnergy_double(new double);
    //std::unique_ptr<double3> h_boxSize_double3(new double3(make_double3(boxInfo.boxSize[0], boxInfo.boxSize[1], boxInfo.boxSize[2])));

    //*totalPEnergy_double = double(totalPEnergy);

    //// Use bulk copy operations and eliminate overhead from individual element access
    //for (size_t i = 0; i < atomPositions.size(); ++i) {
    //    atomPositions_double3[i] = make_double3(atomPositions[i][0], atomPositions[i][1], atomPositions[i][2]);
    //    forces_double3[i] = make_double3(totalForces[i][0], totalForces[i][1], totalForces[i][2]);
    //}

    //// Bulk copy BondParams using std::copy
    //std::copy(bondParams.begin(), bondParams.end(), bondParams_ptr.get());

    //// Launch the kernel with converted data
    //launchKernelBondForces(atomPositions_double3.get(), bondParams_ptr.get(), forces_double3.get(),
    //    totalPEnergy_double.get(), h_boxSize_double3.get(), atomPositions.size(), bondParams.size());

    //// Convert results from double3 back to Coords3D after kernel execution
    //for (size_t i = 0; i < atomPositions.size(); ++i) {
    //    totalForces[i] = Coords3D{ forces_double3[i].x, forces_double3[i].y, forces_double3[i].z };
    //}

    //totalPEnergy = double(*totalPEnergy_double);


    //// Convert double3 results back to Coords3D
    //for (size_t i = 0; i < atomPositions.size(); ++i) {
    //    totalForces[i] = forces_ptr[i];
    //}

    //totalPEnergy = *totalPEnergy_ptr;


    //// Convert Coords3D (host) to double3 (device)
    //Coords3D* atomPositions_ptr = new Coords3D[atomPositions.size()];
    //Coords3D* forces_ptr = new Coords3D[totalForces.size()];
    //BondParams* bondParams_ptr = new BondParams[bondParams.size()];
    //double* totalPEnergy_ptr = new double;
    //Coords3D* h_boxSize_ptr = new Coords3D;

    //*h_boxSize_ptr = boxInfo.boxSize;

    //*totalPEnergy_ptr = totalPEnergy;
    //*atomPositions_ptr = atomPositions[0];
    //*forces_ptr = totalForces[0];
    //*bondParams_ptr = bondParams[0];

    ////for (int i = 0; i < atomPositions.size(); ++i) {
    ////    atomPositions_ptr[i] = atomPositions[i];
    ////    forces_ptr[i] = totalForces[i];
    ////}
    ////for (int j = 0; j < bondParams.size(); ++j) {
    ////    bondParams_ptr[j] = bondParams[j];
    ////}

    // Convert double3 results back to Coords3D
    //for (size_t i = 0; i < atomPositions.size(); ++i) {
    //    totalForces[i] = forces_ptr[i];
    //}

    //totalPEnergy = *totalPEnergy_ptr;

}

void Forces::AddHAngle(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<AngleParams>& angleParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    for (const auto& angle : angleParams) {
        vector<Coords3D> AP = { atomPositions[angle.p1], atomPositions[angle.p2], atomPositions[angle.p3] };
        auto forces = HarmonicAngleForce::calculateForces(AP, angle, totalPEnergy, boxInfo);
        totalForces[angle.p1] += forces[0];
        totalForces[angle.p2] += forces[1];
        totalForces[angle.p3] += forces[2];

        //if (abs(totalForces[angle.p1][1]) > 200) {
        //    cout << "";
        //}
        //if (abs(totalForces[angle.p2][1]) > 200) {
        //    cout << "";
        //}
        //if (abs(totalForces[angle.p3][1]) > 200) {
        //    cout << "";
        //}

    }

}


void Forces::AddNonBondElectroPME(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const NonbondedParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const vector<set<int>>& exclusions) {
    NonbondedForce nonbondedForce;
    auto forces = nonbondedForce.calculateForces(atomPositions, totalPEnergy, params, boxInfo, exclusions);

    for (int atom = 0; atom < atomPositions.size(); atom++) {
        totalForces[atom] += forces[atom];
    }
}

