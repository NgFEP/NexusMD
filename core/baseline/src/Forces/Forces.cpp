#include "Forces.h"
#include "SystemXMLParser.h"
#include "StateXMLParser.h"
#include <utility>

using namespace std;
using namespace BaseLine;


std::vector<Coords3D> Forces::AddPeriodicTorsion(vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<TorsionParameters>& torsionParams) {
    //std::vector<Coords3D> totalForces(atomPositions.size(), Coords3D(0, 0, 0)); // Initialize totalForces with the size of atomPositions

    for (const auto& torsion : torsionParams) {
        //AP is a vector of 4 atoms positions involved in the torsion 
        std::vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
                                     Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
                                     Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
                                     Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };
        auto forces = PeriodicTorsionForce::calculateForces(AP, torsion);
        totalForces[torsion.p1] += forces[0];
        totalForces[torsion.p2] += forces[1];
        totalForces[torsion.p3] += forces[2];
        totalForces[torsion.p4] += forces[3];
    }

    return totalForces;
}

