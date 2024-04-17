#ifndef INITIALIZER_H
#define INITIALIZER_H

#include "Coords3D.h"
#include <vector>
//#include <string>
//#include <utility>//for std::pair

namespace BaseLine {
    class Initializer {
    public:
        Initializer() = default;
        //pair<std::vector<Coords3D>, std::vector<Coords3D>> InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions);
        void InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
    private:
        //std::vector<Coords3D> _totalForces;
        //std::vector<Coords3D> _velocities;
    };

} // namespace BaseLine

#endif // INITIALIZER_H
