#ifndef VERLETINTEGRATION_H
#define VERLETINTEGRATION_H

#include "Coords3D.h"
#include <vector>

namespace BaseLine {

    class VerletIntegration {
    public:
        VerletIntegration();
        ~VerletIntegration();

        static void advance(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& StepNum, double& StepSize);

    private:
        static std::vector<double> inverseMasses;
        static std::vector<Coords3D> UpdatedAtomPositions;
        //void initializeInverseMasses(const std::vector<double>& masses);
        //double getDeltaT() const;
        //int getTimeStep() const;
        //void setTimeStep(int step);
    };

}

#endif // VERLETINTEGRATION_H
