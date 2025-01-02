#include "CudaForceMapper.h"

using namespace std;

CudaForceMapper::CudaForceMapper() {};
CudaForceMapper::~CudaForceMapper() {};

void CudaForceMapper::cudaAllocateBonds(vector<Residues>& _residues, vector<CudaBonds>& _cudaBonds, const vector<vector<int>>& blockResidues) {
    for (const auto& block : blockResidues) { // Iterate over each block
        int cumulativeSize = 0; // Tracks cumulative size of AllAtomsIndices from previous residues

        CudaBonds cudaBond;

        for (int residueIdx : block) { // Process each residue in the block
            const Residues& res = _residues[residueIdx];

            // Copy AllAtomsIndices
            cudaBond.AllAtomsIndices.insert(cudaBond.AllAtomsIndices.end(),
                res.AllAtomsIndices.begin(),
                res.AllAtomsIndices.end());

            // Copy HAtomsIndices
            cudaBond.HAtomsIndices.insert(cudaBond.HAtomsIndices.end(),
                res.HAtomsIndices.begin(),
                res.HAtomsIndices.end());

            // Copy NonHAtomsIndices
            cudaBond.NonHAtomsIndices.insert(cudaBond.NonHAtomsIndices.end(),
                res.NonHAtomsIndices.begin(),
                res.NonHAtomsIndices.end());

            // Copy AllBondsIndices if it exists
            if (res.AllBondsIndices.has_value()) {
                const auto& bonds = res.AllBondsIndices.value();
                for (const auto& bond : bonds) {
                    CudaBondInfo cudaBondInfo = {
                        bond.d,
                        bond.k,
                        bond.bondInx,                        // Bond index           
                        bond.p1InxLocal + cumulativeSize,         // Adjusted p1Inx
                        bond.p2InxLocal + cumulativeSize,         // Adjusted p2Inx
                        bond.p1InxGlobal,
                        bond.p2InxGlobal,
                        (res.resName == "WAT"),
                        bond.p1InRes,                        // Is p1 in residue
                        bond.p2InRes                         // Is p2 in residue

                    };
                    cudaBond.AllBondsIndices.push_back(cudaBondInfo);
                }
            }

            // Copy HBondsIndices if it exists
            if (res.HBondsIndices.has_value()) {
                const auto& hbonds = res.HBondsIndices.value();
                for (const auto& bond : hbonds) {
                    CudaBondInfo cudaBondInfo = {
                        bond.d,
                        bond.k,
                        bond.bondInx,                        // Bond index           
                        bond.p1InxLocal + cumulativeSize,         // Adjusted p1Inx
                        bond.p2InxLocal + cumulativeSize,         // Adjusted p2Inx
                        bond.p1InxGlobal,
                        bond.p2InxGlobal,
                        (res.resName == "WAT"),
                        bond.p1InRes,                        // Is p1 in residue
                        bond.p2InRes                         // Is p2 in residue
                    };
                    cudaBond.HBondsIndices.push_back(cudaBondInfo);
                }
            }

            // Copy NonHBondsIndices if it exists
            if (res.NonHBondsIndices.has_value()) {
                const auto& nonHBonds = res.NonHBondsIndices.value();
                for (const auto& bond : nonHBonds) {
                    CudaBondInfo cudaBondInfo = {
                        bond.d,
                        bond.k,
                        bond.bondInx,                        // Bond index           
                        bond.p1InxLocal + cumulativeSize,         // Adjusted p1Inx
                        bond.p2InxLocal + cumulativeSize,         // Adjusted p2Inx
                        bond.p1InxGlobal,
                        bond.p2InxGlobal,
                        (res.resName == "WAT"),
                        bond.p1InRes,                        // Is p1 in residue
                        bond.p2InRes                         // Is p2 in residue
                    };
                    cudaBond.NonHBondsIndices.push_back(cudaBondInfo);
                }
            }

            // Update cumulativeSize for the next residue
            cumulativeSize += res.AllAtomsIndices.size();
        }

        // Add to _cudaBonds
        _cudaBonds.push_back(move(cudaBond));
    }
}
