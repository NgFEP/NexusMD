#ifndef SEGMENT_FORCE_MAPPER_H
#define SEGMENT_FORCE_MAPPER_H

#include "SystemXMLParser.h"
#include "PDBSegmentParser.h" // Include for PResidues, Connection, WaterMols
#include <vector>
#include <optional>
#include <set>
#include <algorithm> 

// Struct for storing remained bonds
struct RemainedBonds {
    std::set<int> atomIDs; // Unique set of atom indices
    std::optional<std::vector<int>> BondsIndices; // List of all bond indices
};

class SegmentForceMapper {
public:
    // Constructor
    SegmentForceMapper();
    ~SegmentForceMapper();

    void allocateBonds(const std::vector<BondParams>& bondParams, std::vector<PResidues>& residues, std::vector<Connection>& connections, std::vector<WaterMols>& waterMols, RemainedBonds& remainedBonds);


private:
    // Member variables
    int waterMolCounter = 0;
    int connectionCounter = 0;
    int residueCounter = 0;
    // Internal methods
    bool bondInRange(int p1, int p2, int low, int high) const;
    bool bondBetweenAtoms(const BondParams& bond, int atom1, int atom2) const;
    void allocateToResidue(PResidues& residue, const int& bondIndex, const BondParams& bond);
};

#endif // SEGMENT_FORCE_MAPPER_H
