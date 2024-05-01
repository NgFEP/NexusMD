#ifndef HARMONICBONDFORCE_H
#define HARMONICBONDFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //this contains definitions for BondParameters


namespace BaseLine {

    class HarmonicBondForce {
    public:
        static std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, const HBondParams& params, double& totalPEnergy);


    private:

    };

} // namespace BaseLine

#endif // HARMONICBONDFORCE_H