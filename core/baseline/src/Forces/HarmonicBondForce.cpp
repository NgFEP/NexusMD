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


vector<Coords3D> HarmonicBondForce::calculateForces(const vector<Coords3D>& ap, const BondParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    const Coords3D& pos1 = ap[0];
    const Coords3D& pos2 = ap[1];


    // To compute the vector from pos1 to pos2
    Coords3D delta = PeriodicBoundaryCondition::minimumImageVector(pos2 ,pos1, boxInfo);


    // To calculate the distance between the two particles
    double r = sqrt(delta.dot(delta));

    // To compute the energy contribution (this would usually be used elsewhere)
    double deltaIdeal = r - params.d;

    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
    totalPEnergy += energy;

    //std::cout << std::fixed << std::setprecision(16) << "energy: " << energy << std::endl;

    // Compute the derivative of the energy with respect to the distance
    double dEdR = params.k * deltaIdeal;

    // To normalize the delta vector and scale by dEdR
    if (r > 0) {
        delta *= (dEdR / r);
    }
    else {
        delta = Coords3D(0, 0, 0);
    }

    // To calculate the forces on each particle
    Coords3D force1 = delta;       
    Coords3D force2 = -delta;      

    vector<Coords3D> forces = { force1, force2 };
    return forces;
}