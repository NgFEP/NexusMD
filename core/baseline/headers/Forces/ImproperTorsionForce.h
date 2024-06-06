#ifndef HARMONICANGLEFORCE_H
#define HARMONICANGLEFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for AngleParameters
#include "PeriodicBoundaryCondition.h"


namespace BaseLine {

    class HarmonicAngleForce {
    public:
        static double ImproperTorsionForce::calculateImproperTorsionAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const Coords3D& p4);
        static std::vector<Coords3D> ImproperTorsionForce::calculateForces(const std::vector<Coords3D>& atomPositions, double k, double theta0);



    private:

    };

} // namespace BaseLine

#endif // HARMONICANGLEFORCE_H