#ifndef HARMONICANGLEFORCE_H
#define HARMONICANGLEFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for AngleParameters
#include "PeriodicBoundaryCondition.h"


namespace BaseLine {

    class HarmonicAngleForce {
    public:
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const AngleParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        static double calculateAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

    private:

    };

} // namespace BaseLine

#endif // HARMONICANGLEFORCE_H