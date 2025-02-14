#include "Exclusions.h"

using namespace std;
using namespace BaseLine;


void Exclusions::createExclusions(int& numParticles, vector<BondParams>& bonds, vector<set<int>>& exclusions, int& bondCutoff) {
    if (bondCutoff < 1)
        return;

    // Check for illegal particle indices
    for (const auto& bond : bonds) {
        if (bond.p1 < 0 || bond.p2 < 0 || bond.p1 >= numParticles || bond.p2 >= numParticles)
            throw std::runtime_error("createExclusions: Illegal particle index in list of bonds");
    }

    // To initialize exclusions and bonded12 vectors
    exclusions.resize(numParticles);
    vector<set<int>> bonded12(numParticles);

    // Process bonds
    for (const auto& bond : bonds) {
        int p1 = bond.p1;
        int p2 = bond.p2;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }

    // To expand exclusions based on bondCutoff
    for (int level = 0; level < bondCutoff - 1; level++) {
        vector<set<int>> currentExclusions = exclusions;
        for (int i = 0; i < numParticles; i++) {
            for (int j : currentExclusions[i]) {
                for (int bonded : bonded12[i]) {
                    if (bonded != j) {// To ensure that an atom is not added to its own exclusion list
                        exclusions[j].insert(bonded);
                    }
                }
            }
        }
    }
}

