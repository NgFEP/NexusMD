#include "Forces.h"
#include "SystemXMLParser.h"
#include "StateXMLParser.h"
#include <utility>

using namespace std;
using namespace BaseLine;


std::vector<Coords3D> Forces::AddPeriodicTorsion(vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<TorsionParameters>& torsionParams) {
    //std::vector<Coords3D> totalForces(atomPositions.size(), Coords3D(0, 0, 0)); // Initialize totalForces with the size of atomPositions

    //for (const auto& torsion : torsionParams) {
    //    //AP is a vector of 4 atoms positions involved in the torsion 
    //    std::vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
    //                                 Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
    //                                 Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
    //                                 Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };
    //    auto forces = PeriodicTorsionForce::calculateForces(AP, torsion);
    //    totalForces[torsion.p1] += forces[0];
    //    totalForces[torsion.p2] += forces[1];
    //    totalForces[torsion.p3] += forces[2];
    //    totalForces[torsion.p4] += forces[3];
    //}


    //bug found: Forces are calculated from scratch every step and we need to make sure that the vector gets equal to zero on each step before adding up new calculated forces.
    size_t numAtoms = atomPositions.size();
    totalForces.assign(numAtoms, Coords3D(0, 0, 0));


    //For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0]
    //reason for error for (const auto& torsion : torsionParams[0]) { 
    //torsionParams[0] refers to a single instance of TorsionParameters, not a collection that I can iterate over with a range - based for loop.
    if (!torsionParams.empty()) {
        const auto& torsion = torsionParams[0];

        //AP is a vector of 4 atoms positions involved in the torsion 
        std::vector<Coords3D> AP = { Coords3D(atomPositions[torsion.p1][0], atomPositions[torsion.p1][1], atomPositions[torsion.p1][2]),
                                     Coords3D(atomPositions[torsion.p2][0], atomPositions[torsion.p2][1], atomPositions[torsion.p2][2]),
                                     Coords3D(atomPositions[torsion.p3][0], atomPositions[torsion.p3][1], atomPositions[torsion.p3][2]),
                                     Coords3D(atomPositions[torsion.p4][0], atomPositions[torsion.p4][1], atomPositions[torsion.p4][2]) };

        auto forces = PeriodicTorsionForce::calculateForces(AP, torsion);
        // Assuming forcesAndAngle.second is the vector<Coords3D> of forces
        totalForces[torsion.p1] += forces[0];
        totalForces[torsion.p2] += forces[1];
        totalForces[torsion.p3] += forces[2];
        totalForces[torsion.p4] += forces[3];
    }



    return totalForces;
}

