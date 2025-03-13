#include "EnsembleParameters.h"

using namespace std;
using namespace BaseLine;

// Function to calculate DOF for the system based on constraints

// Constructor
EnsembleParameters::EnsembleParameters() {
    // Initialization if needed
}

// Destructor
EnsembleParameters::~EnsembleParameters() {
    // Cleanup if needed
}




int EnsembleParameters::numDOF(int& numAtoms, bool& removeCOM, const vector<Constraint>& constraints) {
    // Start with the unconstrained DOF (3N - 3)
    int dof = 3*numAtoms;  // Total DOF for the unconstrained system

    // Apply constraints
    for (const auto& constraint : constraints) {
        if (constraint.isHBond) {
            dof -= constraint.numHBonds;  // Subtract number of O-H bond constraints
        }
        if (constraint.isHAngle) {
            dof -= constraint.numHAngles;  // Subtract number of H-O-H angle constraints
        }
    }

    if (removeCOM) {
        dof -= 3;  // Remove the center of mass motion (3 degrees of freedom)
    }

    return dof;// -3 * 500;//specific for waterbox test, modify later
}

// ComputeTemperature implementation
double EnsembleParameters::ComputeTemp(int& numAtoms, double& kineticEnergy, bool& removeCOM, const vector<Constraint>& constraints) {
    // Calculate DOF using the constraints
    int dof = numDOF(numAtoms, removeCOM, constraints);

    // Compute temperature using the equipartition theorem: KE = (DOF / 2) * boltz * T
    return (2.0 * kineticEnergy) / (dof * boltz); // 
}


// Compute total mass from atomic masses
double EnsembleParameters::ComputeTotalMass(const vector<double>& Masses) {
    double TotalMass = 0.0;
    for (double mass : Masses) {
        TotalMass += mass;
    }

    // Convert total mass from amu to grams
    TotalMass = TotalMass / AVOGADRO;

    return TotalMass;
}

// Compute density from total mass and volume
double EnsembleParameters::ComputeDensity(double& TotalMass, double& Volume) {
    if (Volume <= 0.0) {
        cerr << "Error: Volume must be greater than zero." << endl;
        return 0.0;
    }

    // Convert volume from am³ to mL

    double Volume_mL = Volume * NM3_TO_ML;

    return TotalMass / Volume_mL;  // Density in g/mL
}

// Compute volume from box size
double EnsembleParameters::ComputeVolume(const Coords3D& boxSize) {
    return boxSize[0] * boxSize[1] * boxSize[2];  // Volume in nm³
}