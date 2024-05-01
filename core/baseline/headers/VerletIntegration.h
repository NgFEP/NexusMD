// VerletIntegration.h
#ifndef VERLET_INTEGRATION_H
#define VERLET_INTEGRATION_H

#include "Coords3D.h"
#include <vector>

namespace BaseLine {

    class VerletIntegration {
    public:
        VerletIntegration();
        ~VerletIntegration();

        // Advances the positions and velocities of atoms in one timestep
        static void InverseMasses(std::vector<double>& masses, std::vector<double>& inverseMasses);
        //static void Advance(const std::vector<double>&dt, std::vector<Coords3D>&atomPositions, std::vector<Coords3D>&velocities, std::vector<Coords3D>&totalForces, std::vector<double>&inverseMasses, std::vector<Coords3D>&posDelta);
        static void Advance(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& inverseMasses, int& StepNum, const std::vector<double>& dt);

        static void SelectVerletStepSize(std::vector<Coords3D>&velocities, std::vector<Coords3D>&totalForces, std::vector<double>&inverseMasses, std::vector<double>&dt, double errorTol, double maxStepSize);


    private:
        static std::vector<Coords3D> UpdatedAtomPositions;


    };

} // namespace BaseLine

#endif // VERLET_INTEGRATION_H
