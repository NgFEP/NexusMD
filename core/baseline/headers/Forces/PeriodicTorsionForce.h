#ifndef PERIODICTORSIONFORCE_H
#define PERIODICTORSIONFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h"
#include "PeriodicBoundaryCondition.h"

namespace BaseLine {


    class PeriodicTorsionForce {
    public:
        PeriodicTorsionForce();
        virtual ~PeriodicTorsionForce();

        // To Calculate the distance between two points
        static double atom_distance(const Coords3D& coords1, const Coords3D& coords2);

        // To Calculate torsion angle (in radians) between four points
        static double torsion_angle(const Coords3D& r1, const Coords3D& r2, const Coords3D& r3, const Coords3D& r4, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        // To Calculate forces for a given torsion
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const PTorsionParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
    };

} // namespace BaseLine

#endif // PERIODICTORSIONFORCE_H
