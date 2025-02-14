#ifndef INITIALIZER_H
#define INITIALIZER_H

#include "Coords3D.h"
#include <vector>


namespace BaseLine {
    class Initializer {
    public:
        Initializer() = default;
        void InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
    private:

    };

} // namespace BaseLine

#endif // INITIALIZER_H
