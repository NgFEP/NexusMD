#ifndef CUDAINITIALIZER_H
#define CUDAINITIALIZER_H

#include "Coords3D.h"
#include <vector>
//#include <string>
//#include <utility>//for std::pair

namespace Cuda {
    class Initializer {
    public:
        Initializer() = default;
        //pair<std::vector<Coords3D>, std::vector<Coords3D>> InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions);
        void InitializeForcesAndVelocities(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
    private:
        //std::vector<Coords3D> _totalForces;
        //std::vector<Coords3D> _velocities;
    };

} // namespace Cuda

#endif // CUDAINITIALIZER_H
