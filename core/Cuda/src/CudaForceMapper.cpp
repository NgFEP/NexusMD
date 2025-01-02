#include "CudaForceMapper.h"

using namespace std;

CudaForceMapper::CudaForceMapper() {};
CudaForceMapper::~CudaForceMapper() {};

void CudaForceMapper::cudaAllocateBonds(vector<Residues>& _residues, vector<CudaBonds>& _cudaBonds, vector<int>& _startResidues, vector<int>& _endResidues) {
    for (int i = 0; i < _startResidues.size(); ++i) {
        int start = _startResidues[i];
        int end = _endResidues[i];


        int cumulativeSize = 0; // Tracks cumulative size of AllAtomsIndices from previous residues

        CudaBonds cudaBond;

        for (int j = start; j <= end; ++j) {
            const Residues& res = _residues[j];

            // Copy AllAtomsIndices
            cudaBond.AllAtomsIndices.insert(cudaBond.AllAtomsIndices.end(),res.AllAtomsIndices.begin(),res.AllAtomsIndices.end());

            // Copy HAtomsIndices
            cudaBond.HAtomsIndices.insert(cudaBond.HAtomsIndices.end(), res.HAtomsIndices.begin(), res.HAtomsIndices.end());

            // Copy NonHAtomsIndices
            cudaBond.NonHAtomsIndices.insert(cudaBond.NonHAtomsIndices.end(), res.NonHAtomsIndices.begin(), res.NonHAtomsIndices.end());

            // Copy AllBondsIndices if it exists
            if (res.AllBondsIndices.has_value()) {
                const auto& bonds = res.AllBondsIndices.value();
                for (const auto& bond : bonds) {
                    CudaBondInfo cudaBondInfo = {
                        bond.bondInx,                        // Bond index           
                        bond.p1Inx + cumulativeSize,         // Adjusted p1Inx
                        bond.p2Inx + cumulativeSize,         // Adjusted p2Inx
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
                        bond.bondInx,
                        bond.p1Inx + cumulativeSize,
                        bond.p2Inx + cumulativeSize,
                        (res.resName == "WAT"),
                        bond.p1InRes,
                        bond.p2InRes
                    };
                    cudaBond.HBondsIndices.push_back(cudaBondInfo);
                }
            }

            // Copy NonHBondsIndices if it exists
            if (res.NonHBondsIndices.has_value()) {
                const auto& nonHBonds = res.NonHBondsIndices.value();
                for (const auto& bond : nonHBonds) {
                    CudaBondInfo cudaBondInfo = {
                        bond.bondInx,
                        bond.p1Inx + cumulativeSize,
                        bond.p2Inx + cumulativeSize,
                        (res.resName == "WAT"),
                        bond.p1InRes,
                        bond.p2InRes
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

