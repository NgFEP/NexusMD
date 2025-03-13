#include "MonteCarloBarostat.h"
#include <cmath>
#include <random>

using namespace std;
using namespace BaseLine;


MonteCarloBarostat::MonteCarloBarostat(double pressure, double temperature, int frequency)
    : _pressure(pressure), _temperature(temperature), _frequency(frequency) {
}


void MonteCarloBarostat::ApplyBarostat(vector<Coords3D>& atomPositions,
    const vector<Molecules>& molecules, int& numMolecules,
    PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // Random number generator for volume scaling
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> randScale(0.95, 1.05); // Small scale factor changes

    double scaleX = randScale(gen);
    double scaleY = randScale(gen);
    double scaleZ = randScale(gen);

    for (int index = 0; index < numMolecules; ++index) {

        // Find the molecular center
        double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
        for (int atom : molecules[index].AllAtomsIndices) {
            centerX += atomPositions[atom][0];
            centerY += atomPositions[atom][1];
            centerZ += atomPositions[atom][2];
        }
        double invNumAtoms = 1.0 / molecules[index].AllAtomsIndices.size();
        centerX *= invNumAtoms;
        centerY *= invNumAtoms;
        centerZ *= invNumAtoms;

        // Compute scaling delta
        double deltaX = centerX * (scaleX - 1.0);
        double deltaY = centerY * (scaleY - 1.0);
        double deltaZ = centerZ * (scaleZ - 1.0);

        // Apply scaling
        for (int atom : molecules[index].AllAtomsIndices) {
            atomPositions[atom][0] += deltaX;
            atomPositions[atom][1] += deltaY;
            atomPositions[atom][2] += deltaZ;
        }
    }

    // Scale the box size
    boxInfo.boxSize[0] *= scaleX;
    boxInfo.boxSize[1] *= scaleY;
    boxInfo.boxSize[2] *= scaleZ;

    //update box boundaries
    PeriodicBoundaryCondition::extractBoxBoundaries(boxInfo);

}
