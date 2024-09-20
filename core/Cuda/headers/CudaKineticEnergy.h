#ifndef CUDAKINETICENERGY_H
#define CUDAKINETICENERGY_H

#include "Coords3D.h"
#include <vector>

namespace Cuda {

    class KineticEnergy {
    public:
        // Calculates kinetic energy for each atom and returns a vector of kinetic energies
        //static void calculateKineticEnergy(const std::vector<Coords3D>& velocities, const std::vector<double>& masses, size_t& numAtoms, std::vector<double>& kineticEnergies, double& totalKEnergy);
        static void calculateKineticEnergy(const std::vector<Coords3D>& velocities, const std::vector<double>& masses, const std::vector<Coords3D>& totalForces, double timeStep, int numAtoms, std::vector<double>& kineticEnergies, double& totalKEnergy);

    };

} // namespace Cuda

#endif // CUDAKINETICENERGY_H
