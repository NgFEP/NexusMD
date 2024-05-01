#include "stdafx.h"
#include "Forces.h"
#include "utility"

using namespace std;
using namespace BaseLine;


void Forces::AddPTorsion(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<PTorsionParams>& torsionParams, double& totalPEnergy) {
    //vector<Coords3D> totalForces(atomPositions.size(), Coords3D(0, 0, 0)); // Initialize totalForces with the size of atomPositions

    for (const auto& torsion : torsionParams) {
        //AP is a vector of 4 atoms positions involved in the torsion 
        vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
                                        Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
                                        Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
                                        Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };
        auto forces = PeriodicTorsionForce::calculateForces(AP, torsion, totalPEnergy);
        totalForces[torsion.p1] += forces[0];
        totalForces[torsion.p2] += forces[1];
        totalForces[torsion.p3] += forces[2];
        totalForces[torsion.p4] += forces[3];
    }



    ////For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0]
    ////reason for error for (const auto& torsion : torsionParams[0]) { 
    ////torsionParams[0] refers to a single instance of TorsionParameters, not a collection that I can iterate over with a range - based for loop.
    //if (!torsionParams.empty()) {
    //    const auto& torsion = torsionParams[0];

    //    //AP is a vector of 4 atoms positions involved in the torsion 
    //    vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
    //                                 Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
    //                                 Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
    //                                 Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };

    //    auto forces = PeriodicTorsionForce::calculateForces(AP, torsion);
    //    // Assuming forcesAndAngle.second is the vector<Coords3D> of forces
    //    totalForces[torsion.p1] += forces[0];
    //    totalForces[torsion.p2] += forces[1];
    //    totalForces[torsion.p3] += forces[2];
    //    totalForces[torsion.p4] += forces[3];
    //}

}


void Forces::AddHBond(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<HBondParams>& bondParams, double& totalPEnergy) {
    for (const auto& bond : bondParams) {
        //if (bond.p1 < atomPositions.size() && bond.p2 < atomPositions.size()) {
        vector<Coords3D> AP = { atomPositions[bond.p1], atomPositions[bond.p2] };
        auto forces = HarmonicBondForce::calculateForces(AP, bond, totalPEnergy);
        totalForces[bond.p1] += forces[0];
        totalForces[bond.p2] += forces[1];
        //}
    }

}

void Forces::AddHAngle(vector<Coords3D>& totalForces, const vector<Coords3D>& atomPositions, const vector<HAngleParams>& angleParams, double& totalPEnergy) {

    for (const auto& angle : angleParams) {
        vector<Coords3D> AP = { atomPositions[angle.p1], atomPositions[angle.p2], atomPositions[angle.p3] };
        auto forces = HarmonicAngleForce::calculateForces(AP, angle, totalPEnergy);
        totalForces[angle.p1] += forces[0];
        totalForces[angle.p2] += forces[1];
        totalForces[angle.p3] += forces[2];
    }

}

