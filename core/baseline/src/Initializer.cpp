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





//pair<vector<Coords3D>, vector<Coords3D>> Initializer::InitializeForcesAndVelocities(const vector<Coords3D>& atomPositions) {
//    size_t numAtoms = atomPositions.size();
//    //vector<Coords3D> totalForces(numAtoms, Coords3D(0, 0, 0));
//    //vector<Coords3D> velocities(numAtoms, Coords3D(0, 0, 0));
//    totalForces.assign(numAtoms, Coords3D(0, 0, 0));
//    velocities.assign(numAtoms, Coords3D(0, 0, 0));
//    // additional initializations here if necessary
//
//    return pair(totalForces, velocities);
//}
void Initializer::InitializeForcesAndVelocities(const vector<Coords3D>& atomPositions, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {
    size_t numAtoms = atomPositions.size();
    totalForces.assign(numAtoms, Coords3D(0, 0, 0));
    velocities.assign(numAtoms, Coords3D(0, 0, 0));

    // Additional initializations here if necessary
}