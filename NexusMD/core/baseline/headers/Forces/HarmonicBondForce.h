#ifndef HARMONICBONDFORCE_H
#define HARMONICBONDFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h"
#include "PeriodicBoundaryCondition.h"


namespace BaseLine {

    class HarmonicBondForce {
    public:
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const BondParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

    private:

    };

} // namespace BaseLine

#endif // HARMONICBONDFORCE_H