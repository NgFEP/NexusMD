#include "Exclusions.h"

using namespace std;
using namespace BaseLine;

//void Exclusions::createExclusions(size_t& numParticles, vector<
// >& bonds,  vector<AngleParams>& angles,  vector<PTorsionParams>& torsions, vector<set<int>>& exclusions) {
//    
//    //if (bondCutoff < 1) return;
//
//
//    //vector<set<int>> bonded12(numParticles);
//
//    // Process bonds
//    for (const auto& bond : bonds) {
//        int p1 = bond.p1;
//        int p2 = bond.p2;
//        exclusions[p1].insert(p2);
//        exclusions[p2].insert(p1);
//        //bonded12[p1].insert(p2);
//        //bonded12[p2].insert(p1);
//    }
//
//    // Process angles
//    for (const auto& angle : angles) {
//        exclusions[angle.p1].insert(angle.p3);
//        exclusions[angle.p3].insert(angle.p1);
//        //bonded12[angle.p1].insert(angle.p3);
//        //bonded12[angle.p3].insert(angle.p1);
//    }
//
//    // Process torsions
//    for (const auto& torsion : torsions) {
//        exclusions[torsion.p1].insert(torsion.p4);
//        exclusions[torsion.p4].insert(torsion.p1);
//        //bonded12[torsion.p1].insert(torsion.p4);
//        //bonded12[torsion.p4].insert(torsion.p1);
//    }
//
//
//}

//void Exclusions::createExclusions(size_t& numParticles, vector<BondParams>& bonds, vector<AngleParams>& angles, vector<PTorsionParams>& torsions, vector<set<int>>& exclusions, int& bondCutoff) {
//    if (bondCutoff < 1) return;
//
//    // Process bonds
//    vector<set<int>> bonded12(numParticles);
//    for (const auto& bond : bonds) {
//        int p1 = bond.p1;
//        int p2 = bond.p2;
//        exclusions[p1].insert(p2);
//        exclusions[p2].insert(p1);
//        bonded12[p1].insert(p2);
//        bonded12[p2].insert(p1);
//    }
//
//    // Process angles
//    for (const auto& angle : angles) {
//        exclusions[angle.p1].insert(angle.p3);
//        exclusions[angle.p3].insert(angle.p1);
//        bonded12[angle.p1].insert(angle.p3);
//        bonded12[angle.p3].insert(angle.p1);
//    }
//
//    // Process torsions
//    for (const auto& torsion : torsions) {
//        exclusions[torsion.p1].insert(torsion.p4);
//        exclusions[torsion.p4].insert(torsion.p1);
//        bonded12[torsion.p1].insert(torsion.p4);
//        bonded12[torsion.p4].insert(torsion.p1);
//    }
//
//    // Expand exclusions based on bondCutoff
//    for (int level = 0; level < bondCutoff - 1; level++) {
//        vector<set<int>> currentExclusions = exclusions;
//        for (size_t i = 0; i < numParticles; i++) {
//            for (int j : currentExclusions[i]) {
//                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
//            }
//        }
//    }
//}



void Exclusions::createExclusions(int& numParticles, vector<BondParams>& bonds, vector<set<int>>& exclusions, int& bondCutoff) {
    if (bondCutoff < 1)
        return;

    // Check for illegal particle indices
    for (const auto& bond : bonds) {
        if (bond.p1 < 0 || bond.p2 < 0 || bond.p1 >= numParticles || bond.p2 >= numParticles)
            throw std::runtime_error("createExclusions: Illegal particle index in list of bonds");
    }

    // Initialize exclusions and bonded12 vectors
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

    // Expand exclusions based on bondCutoff
    for (int level = 0; level < bondCutoff - 1; level++) {
        vector<set<int>> currentExclusions = exclusions;
        for (int i = 0; i < numParticles; i++) {
            for (int j : currentExclusions[i]) {
                for (int bonded : bonded12[i]) {
                    if (bonded != j) {// We need to ensure that we don't add an atom to its own exclusion list
                        exclusions[j].insert(bonded);
                    }
                }
            }
        }
    }
}

