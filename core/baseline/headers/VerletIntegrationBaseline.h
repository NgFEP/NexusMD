#ifndef VERLETINTEGRATION_H
#define VERLETINTEGRATIONBASELINE_H

#include "Coords3D.h"
#include <vector>

namespace BaseLine {

    class VerletIntegrationBaseline {
    public:
        VerletIntegrationBaseline();
        ~VerletIntegrationBaseline();

        static void advance(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& StepNum, double& StepSize);

    private:
        static std::vector<double> inverseMasses;
        static std::vector<Coords3D> UpdatedAtomPositions;
        //void initializeInverseMasses(const vector<double>& masses);
        //double getDeltaT() const;
        //int getTimeStep() const;
        //void setTimeStep(int step);
    };

}

#endif // VERLETINTEGRATIONBASELINE_H
