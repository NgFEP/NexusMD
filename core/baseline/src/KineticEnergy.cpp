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


//// Function to calculate kinetic energy for each atom
//void KineticEnergy::calculateKineticEnergy(const vector<Coords3D>& velocities, const vector<double>& masses, size_t& numAtoms, vector<double>& kineticEnergies, double& totalKEnergy) {
//
//    for (size_t i = 0; i < numAtoms; ++i) {
//        kineticEnergies[i] = 0.5 * masses[i] * velocities[i].dot(velocities[i]);
//        totalKEnergy += kineticEnergies[i];
//    }
//}


// Calculate kinetic energy with time shift
void KineticEnergy::calculateKineticEnergy(const vector<Coords3D>& velocities, const vector<double>& masses, const vector<Coords3D>& totalForces, double timeStep, int numAtoms, vector<double>& kineticEnergies, double& totalKEnergy) {
    const double timeShift = 0.5 * timeStep; // Half of the time step
    //const double scale = timeShift / (double)0x100000000;

    //kineticEnergies.clear();
    //kineticEnergies.resize(numAtoms, 0.0);
    //totalKEnergy = 0.0;

    for (size_t i = 0; i < numAtoms; ++i) {
        Coords3D shiftedVelocity = velocities[i];
        // Apply time shift based on the forces and inverse mass (assuming mass is not zero)
        if (masses[i] > 0) {
            double invMass = 1.0 / masses[i];
            shiftedVelocity += { totalForces[i][0] * invMass* timeShift, totalForces[i][1] * invMass* timeShift, totalForces[i][2] * invMass* timeShift };
        }

        // Compute kinetic energy with shifted velocities
        double ke = 0.5 * masses[i] * shiftedVelocity.dot(shiftedVelocity);

        kineticEnergies[i] = ke;
        totalKEnergy += ke;
        if (ke > 200) {
            cout << "";
        }
    }

}


//int main() {
//    // Example usage:
//    size_t numAtoms = 6;
//    vector<Coords3D> velocities(numAtoms), forces(numAtoms);
//    vector<double> masses(numAtoms, 1.0);  // Assuming all masses are 1 for simplicity
//    vector<double> kineticEnergies;
//    double totalKEnergy = 0.0;
//    double timeStep = 0.001;  // Time step in units of your simulation
//
//    // Initialize velocities and forces as needed
//    KineticEnergy::calculateKineticEnergy(velocities, masses, forces, timeStep, numAtoms, kineticEnergies, totalKEnergy);
//
//    cout << "Total Kinetic Energy: " << totalKEnergy << endl;
//
//    return 0;
//}
