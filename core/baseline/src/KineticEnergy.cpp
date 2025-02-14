//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>  // For std::sqrt

#include "KineticEnergy.h"

using namespace std;
using namespace BaseLine;;


// Function to calculate kinetic energy for each atom
// To calculate kinetic energy with time shift
void KineticEnergy::calculateKineticEnergy(const vector<Coords3D>& velocities, const vector<double>& masses, const vector<Coords3D>& totalForces, double timeStep, int numAtoms, vector<double>& kineticEnergies, double& totalKEnergy) {
    const double timeShift = 0.5 * timeStep; // Half of the time step


    for (size_t i = 0; i < numAtoms; ++i) {
        Coords3D shiftedVelocity = velocities[i];
        // To apply time shift based on the forces and inverse mass (assuming mass is not zero)
        if (masses[i] > 0) {
            double invMass = 1.0 / masses[i];
            shiftedVelocity += { totalForces[i][0] * invMass* timeShift, totalForces[i][1] * invMass* timeShift, totalForces[i][2] * invMass* timeShift };
        }

        // To compute kinetic energy with shifted velocities
        double ke = 0.5 * masses[i] * shiftedVelocity.dot(shiftedVelocity);

        kineticEnergies[i] = ke;
        totalKEnergy += ke;

    }

}
