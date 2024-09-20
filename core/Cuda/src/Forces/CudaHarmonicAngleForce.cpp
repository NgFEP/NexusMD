//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "CudaHarmonicAngleForce.h"
#include <cmath>

using namespace std;
using namespace Cuda;


double HarmonicAngleForce::calculateAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // Coords3D v0 = p2 - p1;
    // Coords3D v1 = p2 - p3;

    Coords3D v0 = PeriodicBoundaryCondition::minimumImageVector(p2, p1, boxInfo);
    Coords3D v1 = PeriodicBoundaryCondition::minimumImageVector(p2, p3, boxInfo);


    Coords3D cp = v0.cross(v1);
    double rp = sqrt(cp.dot(cp));
    rp = max(rp, 1.0e-06);
    double r21 = sqrt(v0.dot(v0));
    double r23 = sqrt(v1.dot(v1));
    double dot = v0.dot(v1);
    double cosine = min(max(dot / (r21 * r23), -1.0), 1.0);
    double theta = acos(cosine);
    return theta;
}

vector<Coords3D> HarmonicAngleForce::calculateForces(const vector<Coords3D>& atomPositions, const HAngleParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    double theta = calculateAngle(atomPositions[0], atomPositions[1], atomPositions[2], boxInfo);
    double deltaIdeal = theta - params.a;
    double energy = 0.5 * params.k * deltaIdeal * deltaIdeal;
    totalPEnergy += energy;

    double dEdAngle = params.k * deltaIdeal;

    // Coords3D v0 = atomPositions[1] - atomPositions[0];
    // Coords3D v1 = atomPositions[1] - atomPositions[2];

    Coords3D v0 = PeriodicBoundaryCondition::minimumImageVector(atomPositions[1], atomPositions[0], boxInfo);
    Coords3D v1 = PeriodicBoundaryCondition::minimumImageVector(atomPositions[1], atomPositions[2], boxInfo);



    Coords3D cp = v0.cross(v1);
    double rp = sqrt(cp.dot(cp));
    rp = max(rp, 1.0e-06);
    double r21 = v0.dot(v0);
    double r23 = v1.dot(v1);

    Coords3D force1 = v0.cross(cp) * (dEdAngle / (r21 * rp));
    Coords3D force3 = cp.cross(v1) * (dEdAngle / (r23 * rp));
    Coords3D force2 = -force1 - force3;

    vector<Coords3D> forces = { force1, force2, force3 };
    return forces;
}

