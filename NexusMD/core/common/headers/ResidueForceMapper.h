#ifndef RESIDUE_FORCE_MAPPER_H
#define RESIDUE_FORCE_MAPPER_H

#include "SystemXMLParser.h"
#include "PDBResidueParser.h"
#include <vector>
#include <optional>
#include <set>
#include <algorithm> 
#include <unordered_set>
#include <cuda_runtime.h>

// Struct for storing remained bonds
struct RemainedBonds {
    std::set<int> atomIDs; // Unique set of atom indices
    std::optional<std::vector<int>> BondsIndices; // List of all bond indices
};

class ResidueForceMapper {
public:
    // Constructor
    ResidueForceMapper();
    ~ResidueForceMapper();

    void allocateBonds(const std::vector<BondParams>& bondParams, std::vector<Residues>& residues, RemainedBonds& remainedBonds, int& totalResiduesSize, int& totalBondsInResidues, std::vector<int>& startResidues, std::vector<int>& endResidues);


private:
    // Member variables
    int bondAtomInResidue = 0;
    int nonResAtomInx = -1;
    int isAssigned = 0;// checks bond assignment 2 means bond assigned in completed, 1 partially done, 0 not assigned to any atoms
    std::vector<int> ResAtomIndicesMaxElement;


    // Internal methods
    int calculateResidueMemory(const Residues& residue);

    void bondInResidue(int p1, int p2, const std::vector<int>& allAtomsIndices);
    void allocateToResidue(Residues& residue, const int& bondIndex, const BondParams& bond);
};

#endif // RESIDUE_FORCE_MAPPER_H
