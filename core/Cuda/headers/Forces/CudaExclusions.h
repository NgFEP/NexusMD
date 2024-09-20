#ifndef CUDAEXCLUSIONS_H
#define CUDAEXCLUSIONS_H

#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include "SystemXMLParser.h"

namespace Cuda {

    class Exclusions {
    public:
        //static void createExclusions(std::size_t& numParticles, std::vector<HBondParams>& bonds, std::vector<HAngleParams>& angles, std::vector<PTorsionParams>& torsions, std::vector<std::set<int>>& exclusions, int& bondCutoff);
        static void createExclusions(int& numParticles, std::vector<HBondParams>& bonds, std::vector<std::set<int>>& exclusions, int& bondCutoff);

    private:

    };
}
#endif // CUDAEXCLUSIONS_H
