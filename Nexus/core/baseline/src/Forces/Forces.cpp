//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "Forces.h"
#include "utility"

using namespace std;
using namespace BaseLine;


void Forces::AddPTorsion(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<PTorsionParams>& torsionParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    for (const auto& torsion : torsionParams) {
        //AP is a vector of 4 atoms positions involved in the torsion 
        vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
                                        Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
                                        Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
                                        Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };
        auto forces = PeriodicTorsionForce::calculateForces(AP, torsion, totalPEnergy, boxInfo);
        totalForces[torsion.p1] += forces[0];
        totalForces[torsion.p2] += forces[1];
        totalForces[torsion.p3] += forces[2];
        totalForces[torsion.p4] += forces[3];






    }

}


void Forces::AddHBond(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<BondParams>& bondParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    for (const auto& bond : bondParams) {
        //if (bond.p1 < atomPositions.size() && bond.p2 < atomPositions.size()) {
        vector<Coords3D> AP = { atomPositions[bond.p1], atomPositions[bond.p2] };
        auto forces = HarmonicBondForce::calculateForces(AP, bond, totalPEnergy, boxInfo);
        totalForces[bond.p1] += forces[0];
        totalForces[bond.p2] += forces[1];

    }

}

void Forces::AddHAngle(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<AngleParams>& angleParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    for (const auto& angle : angleParams) {
        vector<Coords3D> AP = { atomPositions[angle.p1], atomPositions[angle.p2], atomPositions[angle.p3] };
        auto forces = HarmonicAngleForce::calculateForces(AP, angle, totalPEnergy, boxInfo);
        totalForces[angle.p1] += forces[0];
        totalForces[angle.p2] += forces[1];
        totalForces[angle.p3] += forces[2];


    }

}


void Forces::AddNonBondElectroPME(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const NonbondedParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const vector<set<int>>& exclusions) {
    NonbondedForce nonbondedForce;
    auto forces = nonbondedForce.calculateForces(atomPositions, totalPEnergy, params, boxInfo, exclusions);

    for (int atom = 0; atom < atomPositions.size(); atom++) {
        totalForces[atom] += forces[atom];
    }
}

