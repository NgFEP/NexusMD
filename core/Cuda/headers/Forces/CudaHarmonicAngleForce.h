#ifndef CUDAHARMONICANGLEFORCE_H
#define CUDAHARMONICANGLEFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for AngleParameters
#include "PeriodicBoundaryCondition.h"


namespace Cuda {

    class HarmonicAngleForce {
    public:
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const HAngleParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        static double calculateAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

    private:

    };

} // namespace Cuda

#endif // CUDAHARMONICANGLEFORCE_H