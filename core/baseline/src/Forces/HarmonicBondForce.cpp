//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "HarmonicBondForce.h"
#include <cmath>
#include <iomanip>//to print with 16 decimal  

using namespace std;
using namespace BaseLine;


vector<Coords3D> HarmonicBondForce::calculateForces(const vector<Coords3D>& ap, const HBondParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // assert(params.p1 < atomPositions.size() && params.p2 < atomPositions.size());

    //bug: this step of assigning positions of the specific potential has been already taken care of in the forces step 
    //const Coords3D& pos1 = atomPositions[params.p1];
    //const Coords3D& pos2 = atomPositions[params.p2];

    const Coords3D& pos1 = ap[0];
    const Coords3D& pos2 = ap[1];


    // Compute the vector from pos1 to pos2
    // Coords3D delta = pos2 - pos1;
    Coords3D delta = PeriodicBoundaryCondition::minimumImageVector(pos2 ,pos1, boxInfo);

    //std::cout << std::fixed << std::setprecision(16) << "delta: " << delta << std::endl;

    // Calculate the distance between the two particles
    double r = sqrt(delta.dot(delta));
    //std::cout << std::fixed << std::setprecision(16) << "delta[0]: " << delta[0] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "delta[0]*delta[0]: " << delta[0] * delta[0] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "delta[1]: " << delta[1] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "delta[1]*delta[1]: " << delta[1] * delta[1] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "delta[2]: " << delta[2] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "delta[2]*delta[2]: " << delta[2] * delta[2] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "r2: " << delta[0] * delta[0]+ delta[1] * delta[1]+ delta[2] * delta[2] << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "r: " << r << std::endl;

    // Compute the energy contribution (this would usually be used elsewhere)
    double deltaIdeal = r - params.d;
    //std::cout << std::fixed << std::setprecision(16) << "bondParams distance: " << params.d << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "deltaIdeal: " << deltaIdeal << std::endl;
    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
    totalPEnergy += energy;

    //std::cout << std::fixed << std::setprecision(16) << "energy: " << energy << std::endl;

    // Compute the derivative of the energy with respect to the distance
    double dEdR = params.k * deltaIdeal;
    //std::cout << std::fixed << std::setprecision(16) << "bondParams force constant: " << params.k << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "dEdR: " << dEdR << std::endl;

    // Normalize the delta vector and scale by dEdR
    if (r > 0) {
        delta *= (dEdR / r);
    }
    else {
        delta = Coords3D(0, 0, 0);  // Prevent division by zero
    }
    //std::cout << std::fixed << std::setprecision(16) << "normalized delta: " << delta << std::endl;

    // Calculate the forces on each particle
    Coords3D force1 = delta;       // Force on p1
    Coords3D force2 = -delta;      // Force on p2, equal and opposite

    vector<Coords3D> forces = { force1, force2 };
    return forces;
}