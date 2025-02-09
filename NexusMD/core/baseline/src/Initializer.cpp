//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "Initializer.h"
#include "utility"

using namespace std;
using namespace BaseLine;


void Initializer::InitializeForcesAndVelocities(const vector<Coords3D>& atomPositions, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {
    int numAtoms = atomPositions.size();
    totalForces.assign(numAtoms, Coords3D(0, 0, 0));
    velocities.assign(numAtoms, Coords3D(0, 0, 0));

    // Additional initializations here if necessary
}