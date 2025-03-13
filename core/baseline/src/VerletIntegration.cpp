//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "Coords3D.h"
#include <cstring>
#include <cstdio>
#include "StateXMLParser.h"
#include "VerletIntegration.h"
#include "PeriodicBoundaryCondition.h"
#include <cmath> 
#include <algorithm>


using namespace BaseLine;
using namespace std;



// To define the static member variable with namespace qualification
vector<Coords3D> VerletIntegration::UpdatedAtomPositions;
VerletIntegration::VerletIntegration() {
    // Constructor implementation
}

VerletIntegration::~VerletIntegration() {
    // Destructor implementation
}



//modified version Nexus modified with cuda openmm method
void VerletIntegration::Advance(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& inverseMasses, int& StepNum, const vector<double>& dt, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    double dtPos = dt[1];
    double dtVel = 0.5 * (dt[0] + dt[1]);
    double scale = dtVel;// / (double)0x100000000;
    
    int numberOfAtoms = atomPositions.size();

    UpdatedAtomPositions.assign(numberOfAtoms, Coords3D(0, 0, 0));

    // To perform the integration
    for (int i = 0; i < numberOfAtoms; ++i) {
        if (inverseMasses[i] != 0.0) {
            for (int j = 0; j < 3; ++j) {
                //note that the following calculated velocity is at t+StepSize/2 and not t+StepSize
                velocities[i][j] += inverseMasses[i] * totalForces[i][j] * scale;


                //note that the following calculated atompositions is at t+StepSize
                UpdatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j] * dtPos;

                // For PBC: Wrap coordinates to stay within box dimensions
                if (UpdatedAtomPositions[i][j] < boxInfo.lb[j]) {
                    UpdatedAtomPositions[i][j] += boxInfo.boxSize[j];
                }
                else if (UpdatedAtomPositions[i][j] > boxInfo.ub[j]) {
                    UpdatedAtomPositions[i][j] -= boxInfo.boxSize[j];
                }
            }
        }
    }
    double oneOverDt = 1.0 / dt[1];

    for (int i = 0; i < numberOfAtoms; ++i) {
        if (inverseMasses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {

                atomPositions[i][j] = UpdatedAtomPositions[i][j];
            }
    }

}

// dynamically adjust the integration step size based on the forces experienced by the particles and their masses
void VerletIntegration::SelectVerletStepSize(vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& inverseMasses, vector<double>& dt, double errorTol, double maxStepSize) {
    double totalError = 0.0;
    double scale = 1.0; // / (double)0x100000000;

    for (int i = 0; i < velocities.size(); ++i) {
        double invMass = inverseMasses[i];
        if (invMass != 0.0) {
            double err = (scale * totalForces[i][0]) * (scale * totalForces[i][0]) +
                (scale * totalForces[i][1]) * (scale * totalForces[i][1]) +
                (scale * totalForces[i][2]) * (scale * totalForces[i][2]);
            err *= invMass * invMass;
            totalError += err;
        }
    }
    
    totalError = sqrt(totalError / (3 * velocities.size()));
    double newStepSize = sqrt(errorTol / totalError);
    double oldStepSize = dt[1];
    if (oldStepSize > 0) {
        newStepSize = min(newStepSize, oldStepSize * 2.0);
    }
    if (newStepSize > oldStepSize && newStepSize < 1.1 * oldStepSize) {
        newStepSize = oldStepSize;
    }
    if (newStepSize > maxStepSize) {
        newStepSize = maxStepSize;
    }
    dt[1] = newStepSize;
}