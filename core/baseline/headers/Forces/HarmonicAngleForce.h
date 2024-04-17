#ifndef HARMONICANGLEFORCE_H
#define HARMONICANGLEFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for AngleParameters


namespace BaseLine {

    class HarmonicAngleForce {
    public:
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const HAngleParams& params);
        static double calculateAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3);

    private:

    };

} // namespace BaseLine

#endif // HARMONICANGLEFORCE_H