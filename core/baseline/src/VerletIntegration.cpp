#include "Coords3D.h"
#include <vector>
#include <cstring>
#include <sstream>
#include <cstdio>
#include "StateXMLParser.h"
#include "VerletIntegration.h"

using namespace BaseLine;
using namespace std;

#include "VerletIntegration.h"
#include <cmath> // Include for mathematical functions like sqrt()

// Define the static member variable with namespace qualification
std::vector<double> VerletIntegration::inverseMasses;

VerletIntegration::VerletIntegration() {
    // Constructor implementation
}

VerletIntegration::~VerletIntegration() {
    // Destructor implementation
}

void VerletIntegration::advance(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& StepNum, double& StepSize) {
    int numberOfAtoms = atomPositions.size();
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
                velocities[i][j] += inverseMasses[i] * totalForces[i][j] * StepSize;
                atomPositions[i][j] += velocities[i][j] * StepSize;
            }
        }
    }
}


