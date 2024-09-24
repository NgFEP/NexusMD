#ifndef HARMONICANGLEFORCE_H
#define HARMONICANGLEFORCE_H

#include "CudaCoords3D.h"
#include <vector>
#include "CudaSystemXMLParser.h" //this contains definitions for AngleParameters
#include "CudaPeriodicBoundaryCondition.h"


namespace Cuda {

    class HarmonicAngleForce {
    public:
        static double ImproperTorsionForce::calculateImproperTorsionAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const Coords3D& p4);
        static std::vector<Coords3D> ImproperTorsionForce::calculateForces(const std::vector<Coords3D>& atomPositions, double k, double theta0);



    private:

    };

} // namespace Cuda

#endif // HARMONICANGLEFORCE_H