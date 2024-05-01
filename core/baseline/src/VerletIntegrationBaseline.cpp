#include "stdafx.h"
#include "Coords3D.h"
#include <cstring>
#include <cstdio>
#include "StateXMLParser.h"
#include "VerletIntegrationBaseline.h"
#include <cmath> // Include for mathematical functions like sqrt()

using namespace BaseLine;
using namespace std;



// Define the static member variable with namespace qualification
vector<double> VerletIntegrationBaseline::inverseMasses;
vector<Coords3D> VerletIntegrationBaseline::UpdatedAtomPositions;
VerletIntegrationBaseline::VerletIntegrationBaseline() {
    // Constructor implementation
}

VerletIntegrationBaseline::~VerletIntegrationBaseline() {
    // Destructor implementation
}

void VerletIntegrationBaseline::advance(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& StepNum, double& StepSize) {
    int numberOfAtoms = atomPositions.size();

    UpdatedAtomPositions.assign(numberOfAtoms, Coords3D(0, 0, 0));

    // Initialize inverse masses only once to prevent overcomputation of inverseMasses for every steps
    if (StepNum == 0) {
        inverseMasses.resize(masses.size());
        for (size_t i = 0; i < masses.size(); ++i) {
            inverseMasses[i] = (masses[i] == 0.0) ? 0.0 : 1.0 / masses[i];
        }
    }

    // Perform the integration
    for (int i = 0; i < numberOfAtoms; ++i) {
        if (masses[i] != 0.0) {
            for (int j = 0; j < 3; ++j) {
                //note that the following calculated velocity is at t+StepSize/2 and not t+StepSize
                velocities[i][j] += inverseMasses[i] * totalForces[i][j] * StepSize/2.0;
                //note that the following calculated atompositions is at t+StepSize
                UpdatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j] * StepSize;
            }
        }
    }
    for (int i = 0; i < numberOfAtoms; ++i) {
        if (masses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {
                //now these velocities and atomPositions are at the next step of t+StepSize 
                velocities[i][j] = (UpdatedAtomPositions[i][j] - atomPositions[i][j])/ StepSize;
                atomPositions[i][j] = UpdatedAtomPositions[i][j];
            }
    }
}
