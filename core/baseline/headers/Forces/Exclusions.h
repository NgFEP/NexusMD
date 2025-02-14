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
                static void createExclusions(int& numParticles, std::vector<BondParams>& bonds, std::vector<std::set<int>>& exclusions, int& bondCutoff);

    private:

    };
}
#endif // EXCLUSIONS_H
