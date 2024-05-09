//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "KineticEnergy.h"

using namespace std;
using namespace BaseLine;;



// Function to calculate kinetic energy for each atom
void KineticEnergy::calculateKineticEnergy(const vector<Coords3D>& velocities, const vector<double>& masses, size_t& numAtoms, vector<double>& kineticEnergies, double& totalKEnergy) {

    for (size_t i = 0; i < numAtoms; ++i) {
        kineticEnergies[i] = 0.5 * masses[i] * velocities[i].dot(velocities[i]);
        totalKEnergy += kineticEnergies[i];
    }
}


//int main() {
//    // Example velocities and masses
//    vector<Coords3D> velocities = {
//        Coords3D(2.0, 3.0, 4.0),
//        Coords3D(1.5, 2.5, 3.5)
//    };
//    vector<double> masses = { 2.0, 3.0 };
//
//    // Calculate kinetic energies
//    vector<double> kineticEnergies = KineticEnergy::calculate(velocities, masses);
//
//    // Output the kinetic energies
//    for (double ke : kineticEnergies) {
//        cout << "Kinetic Energy: " << ke << " Joules" << endl;
//    }
//
//    return 0;
//}
