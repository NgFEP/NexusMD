#ifndef CUDAPERIODICTORSIONFORCE_H
#define CUDAPERIODICTORSIONFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for TorsionParameters
#include "PeriodicBoundaryCondition.h"

namespace Cuda {


    class PeriodicTorsionForce {
    public:
        PeriodicTorsionForce();
        virtual ~PeriodicTorsionForce();

        // Calculate the distance between two points
        static double atom_distance(const Coords3D& coords1, const Coords3D& coords2);

        // Calculate torsion angle (in radians) between four points
        static double torsion_angle(const Coords3D& r1, const Coords3D& r2, const Coords3D& r3, const Coords3D& r4, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        // Calculate forces for a given torsion
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const PTorsionParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
    };

} // namespace Cuda

#endif // CUDAPERIODICTORSIONFORCE_H
