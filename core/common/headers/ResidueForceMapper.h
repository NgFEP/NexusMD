#ifndef RESIDUE_FORCE_MAPPER_H
#define RESIDUE_FORCE_MAPPER_H

#include "SystemXMLParser.h"
#include "PDBResidueParser.h" // Include for PResidues, Connection, wResidues
#include <vector>
#include <optional>
#include <set>
#include <algorithm> 

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

    void allocateBonds(const std::vector<BondParams>& bondParams, std::vector<PResidues>& residues, std::vector<WResidues>& wResidues, RemainedBonds& remainedBonds);


private:
    // Member variables
    int wResidueCounter = 0;
    //int connectionCounter = 0;
    //int pResidueCounter = 0;
    int bondAtomInRangeCounter = 0;
    int nonResAtomInx = -1;
    int isAssigned = 0;// checks bond assignment 2 means bond assigned in completed, 1 partially done, 0 not assigned to any atoms

    // Internal methods
    void bondInRangePResidue(int p1, int p2, int low, int high);
    bool bondInRangeWResidue(int p1, int p2, int low, int high);
    bool bondBetweenAtoms(const BondParams& bond, int atom1, int atom2) const;
    void allocateToResidue(PResidues& residue, const int& bondIndex, const BondParams& bond);
};

#endif // RESIDUE_FORCE_MAPPER_H
