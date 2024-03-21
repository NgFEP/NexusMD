#ifndef INITIALIZER_H
#define INITIALIZER_H

#include "Coords3D.h"
#include <vector>
#include <string>
#include <utility>//for std::pair
using namespace std;

namespace BaseLine {
    class Initializer {
    public:
        Initializer() = default;
        //pair<vector<Coords3D>, vector<Coords3D>> InitializeForcesAndVelocities(const vector<Coords3D>& atomPositions);
        void InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
    private:
        vector<Coords3D> totalForces;
        vector<Coords3D> velocities;
    };

} // namespace BaseLine

#endif // INITIALIZER_H
