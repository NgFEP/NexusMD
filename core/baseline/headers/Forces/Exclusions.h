#ifndef EXCLUSIONS_H
#define EXCLUSIONS_H

#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include "SystemXMLParser.h"

namespace BaseLine {

    class Exclusions {
    public:
        //static void createExclusions(std::size_t& numParticles, std::vector<
        // >& bonds, std::vector<AngleParams>& angles, std::vector<PTorsionParams>& torsions, std::vector<std::set<int>>& exclusions, int& bondCutoff);
        static void createExclusions(int& numParticles, std::vector<BondParams>& bonds, std::vector<std::set<int>>& exclusions, int& bondCutoff);

    private:

    };
}
#endif // EXCLUSIONS_H
