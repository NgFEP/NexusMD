#include <cstring>
#include <sstream>
#include <cstdio>

#include "VerletIntegration.h"


void VerletIntegration::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
    vector<Vec3>& velocities,
    vector<Vec3>& forces, vector<double>& masses, double tolerance) {
    // first-time-through initialization

    int numberOfAtoms = system.getNumParticles();
    if (getTimeStep() == 0) {
        // invert masses

        for (int ii = 0; ii < numberOfAtoms; ii++) {
            if (masses[ii] == 0.0)
                inverseMasses[ii] = 0.0;
            else
                inverseMasses[ii] = 1.0 / masses[ii];
        }
    }

    // Perform the integration.

    for (int i = 0; i < numberOfAtoms; ++i) {
        if (masses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {
                velocities[i][j] += inverseMasses[i] * forces[i][j] * getDeltaT();
                xPrime[i][j] = atomCoordinates[i][j] + velocities[i][j] * getDeltaT();
            }
    }
}
