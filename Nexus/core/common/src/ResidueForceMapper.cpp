#include "ResidueForceMapper.h"

using namespace std;

ResidueForceMapper::ResidueForceMapper() {};
ResidueForceMapper::~ResidueForceMapper() {};


int ResidueForceMapper::calculateResidueMemory(const Residues& residue) {
    // Size of mandatory vector members
    int atomMemory = residue.AllAtomsIndices.size() * sizeof(int);
    int hAtomMemory = residue.HAtomsIndices.size() * sizeof(int);
    int nonHAtomMemory = residue.NonHAtomsIndices.size() * sizeof(int);

    // Size of optional vector members (check if they have a value before accessing)
    int bondMemory = residue.AllBondsIndices ? residue.AllBondsIndices->size() * sizeof(int) : 0;
    int hBondMemory = residue.HBondsIndices ? residue.HBondsIndices->size() * sizeof(int) : 0;
    int nonHBondMemory = residue.NonHBondsIndices ? residue.NonHBondsIndices->size() * sizeof(int) : 0;
    int nonResBondAtomMemory = residue.NonResBondAtoms ? residue.NonResBondAtoms->size() * sizeof(int) : 0;

    int angleMemory = residue.AngleIndices ? residue.AngleIndices->size() * sizeof(int) : 0;
    int torsionMemory = residue.TorsionIndices ? residue.TorsionIndices->size() * sizeof(int) : 0;

    // Memory for the string (capacity, since it allocates memory dynamically)
    int staticMemory = residue.resName.capacity() * sizeof(char);

    // Total memory for this residue
    int totalMemory = atomMemory + hAtomMemory + nonHAtomMemory +
        bondMemory + hBondMemory + nonHBondMemory +
        nonResBondAtomMemory + staticMemory;// +angleMemory + torsionMemory;

    return totalMemory;
}


void ResidueForceMapper::allocateBonds(const vector<BondParams>& bondParams,vector<Residues>& residues,RemainedBonds& remainedBonds, int& totalResiduesSize, int& totalBondsInResidues, vector<int>& startResidues, vector<int>& endResidues) {

    set<int> assignedwResidues;
    int numResidues = residues.size();



    for (size_t bondIndex = 0; bondIndex < bondParams.size(); ++bondIndex) {
        const auto& bond = bondParams[bondIndex];
        isAssigned = 0;
        // To check if the bond belongs to any residue
        //if (wResidueCounter==0) {// as soon as water molecules start residue checking ends
        for (size_t i = 0; i < residues.size(); ++i) {
            if (i >= ResAtomIndicesMaxElement.size()) {
                ResAtomIndicesMaxElement.push_back(*std::max_element(residues[i].AllAtomsIndices.begin(), residues[i].AllAtomsIndices.end()));
            }
            if (bond.p1 <= ResAtomIndicesMaxElement[i] || bond.p2 <= ResAtomIndicesMaxElement[i]) {// || ensures connecting bonds are covered in both of the involving residues
                bondInResidue(bond.p1, bond.p2, residues[i].AllAtomsIndices);
                if (bondAtomInResidue == 2) {
                    allocateToResidue(residues[i], bondIndex, bond);
                }
                else if (bondAtomInResidue == 1) {
                    if (!residues[i].NonResBondAtoms) {
                        residues[i].NonResBondAtoms.emplace();
                    }
                    residues[i].NonResBondAtoms->push_back(nonResAtomInx);
                    allocateToResidue(residues[i], bondIndex, bond);

                }
                if (isAssigned == 2) {
                    break;
                }
            }    
        }
        //}
        if (isAssigned==2) {
            continue; // To move to the next bond
        }

        // If bond is not assigned, add it to RemainedBonds
        if (isAssigned==0) {
            remainedBonds.atomIDs.insert(bond.p1);
            remainedBonds.atomIDs.insert(bond.p2);
            if (!remainedBonds.BondsIndices) {
                remainedBonds.BondsIndices.emplace();
            }
            remainedBonds.BondsIndices->push_back(bondIndex);
        }
    }

    // Sorting residues in Descending Order
    std::sort(residues.begin(), residues.end(), [](const Residues& a, const Residues& b) {
        auto sizeA = a.AllBondsIndices ? a.AllBondsIndices->size() : 0;
        auto sizeB = b.AllBondsIndices ? b.AllBondsIndices->size() : 0;
        return sizeA > sizeB; // Change comparison to greater-than
        });



    for (size_t i = 0; i < residues.size(); ++i) {
        int residueMemory = calculateResidueMemory(residues[i]); // To calculate memory for this residue
        residues[i].resMemSize = residueMemory;
        totalResiduesSize += residueMemory;                // To accumulate total memory
    }


    for (int i = 0; i < numResidues; ++i) {
        totalBondsInResidues += residues[i].AllBondsIndices ? residues[i].AllBondsIndices->size() : 0;
    }

    int deviceId;
    cudaGetDevice(&deviceId);

    int sharedMemPerBlock;
    cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);


    // To calculate sizes of static components
    const int sizeOfPositions = sizeof(double3);
    const int sizeOfBondParams = sizeof(BondParams);

    int totalMemoryAllBonds = totalBondsInResidues * (2 * sizeOfPositions + sizeOfBondParams) + totalResiduesSize;
    int numBlocks = 1.1 * totalMemoryAllBonds / sharedMemPerBlock; // 10% margin

    // Example of precomputing residue ranges on the host


    int residueIdx = 0;
    while(residueIdx < numResidues){
        startResidues.push_back(residueIdx);

        int sharedMemoryUsed = 0;
        while (residueIdx < numResidues && (sharedMemoryUsed + residues[residueIdx].AllBondsIndices->size()*(2 * sizeOfPositions + sizeOfBondParams) + residues[residueIdx].resMemSize <= 0.95 * sharedMemPerBlock)) {//0.95 to keep a safe margin of unused sharedMemPerBlock
            sharedMemoryUsed += residues[residueIdx].AllBondsIndices->size() * (2 * sizeOfPositions + sizeOfBondParams) + residues[residueIdx].resMemSize;
            residueIdx++;
        }
        endResidues.push_back(residueIdx-1);
    }
}


// Utility function to check if a bond is within a specific range of PResidue

void ResidueForceMapper::bondInResidue(int p1, int p2, const vector<int>& allAtomsIndices) {
    // To create a hash set for efficient lookup
    std::unordered_set<int> atomSet(allAtomsIndices.begin(), allAtomsIndices.end());

    bool p1InResidue = atomSet.find(p1) != atomSet.end();
    bool p2InResidue = atomSet.find(p2) != atomSet.end();

    if (p1InResidue && p2InResidue) {
        bondAtomInResidue = 2;
        isAssigned = 2;
    }
    else if (p1InResidue) {
        bondAtomInResidue = 1;
        nonResAtomInx = p2;
        isAssigned++;
    }
    else if (p2InResidue) {
        bondAtomInResidue = 1;
        nonResAtomInx = p1;
        isAssigned++;
    }
    else {
        bondAtomInResidue = 0;
    }
}


// Utility function to allocate bonds to a residue
void ResidueForceMapper::allocateToResidue(Residues& residue,const int& bondIndex, const BondParams& bond) {
    // Check if bond.p1 or bond.p2 is a hydrogen atom in this residue
    bool p1IsHydrogen = find(residue.HAtomsIndices.begin(), residue.HAtomsIndices.end(), bond.p1) != residue.HAtomsIndices.end();
    bool p2IsHydrogen = find(residue.HAtomsIndices.begin(), residue.HAtomsIndices.end(), bond.p2) != residue.HAtomsIndices.end();


    // To ensure AllBondsIndices is initialized
    if (!residue.AllBondsIndices) {
        //residue.AllBondsIndices = {};
        residue.AllBondsIndices.emplace();

    }
    residue.AllBondsIndices->push_back(bondIndex);

    // If either atom is a hydrogen atom, classify as a hydrogen bond
    if (p1IsHydrogen || p2IsHydrogen) {
        if (!residue.HBondsIndices) {
            residue.HBondsIndices.emplace();
        }
        residue.HBondsIndices->push_back(bondIndex);
    }
    else { // Otherwise, classify as a non-hydrogen bond
        if (!residue.NonHBondsIndices) {
            residue.NonHBondsIndices.emplace();
        }
        residue.NonHBondsIndices->push_back(bondIndex);
    }
}
