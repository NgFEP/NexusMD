#include "SegmentForceMapper.h"

using namespace std;

SegmentForceMapper::SegmentForceMapper() {};
SegmentForceMapper::~SegmentForceMapper() {};


void SegmentForceMapper::allocateBonds(const vector<BondParams>& bondParams,vector<PResidues>& residues,vector<Connection>& connections,vector<WaterMols>& waterMols,RemainedBonds& remainedBonds) {

    set<int> assignedConnections;
    set<int> assignedWaterMols;

    waterMolCounter = 0;
    connectionCounter = 0;
    residueCounter = 0;

    for (size_t bondIndex = 0; bondIndex < bondParams.size(); ++bondIndex) {
        const auto& bond = bondParams[bondIndex];
        bool isAssigned = false;

        // Check if the bond belongs to any residue
        if (waterMolCounter==0) {// as soon as water molecules start residue checking ends
            for (size_t i = residueCounter; i < residues.size(); ++i) {
                if (bondInRange(bond.p1, bond.p2, residues[i].lowBound, residues[i].highBound)) {
                    allocateToResidue(residues[i], bondIndex, bond);
                    isAssigned = true;
                    break;
                }
            }
        }
        if (isAssigned) {
            continue; // Move to the next bond
        }
        // Check if the bond belongs to any connection
        if (!isAssigned && connectionCounter < connections.size()) {
            for (size_t i = connectionCounter; i < connections.size(); ++i) {
                if (assignedConnections.find(i) == assignedConnections.end() &&
                    bondBetweenAtoms(bond, connections[i].atomC, connections[i].atomN)) {
                    // Store the bond index, not the atom indices
                    connections[i].BondIndex = bondIndex;
                    assignedConnections.insert(i);
                    connectionCounter++;
                    residueCounter++;
                    isAssigned = true;
                    break;
                }
            }
        }
        if (isAssigned) {
            continue; // Move to the next bond
        }
        // Check if the bond belongs to any water molecule
        if (!isAssigned && waterMolCounter < waterMols.size()) {
            for (size_t i = waterMolCounter; i < waterMols.size(); ++i) {
                if (assignedWaterMols.find(i) == assignedWaterMols.end() &&
                    bondInRange(bond.p1, bond.p2, waterMols[i].lowBound, waterMols[i].highBound)) {

                    if (!waterMols[i].BondsIndices) {
                        waterMols[i].BondsIndices.emplace();
                    }

                    // Store the bond index, not the atom indices
                    waterMols[i].BondsIndices->push_back(bondIndex);
                    if (waterMols[i].BondsIndices->size()==2) {//2 bonds alocated to a water molecule
                        assignedWaterMols.insert(i);
                        waterMolCounter++;
                    }
                    isAssigned = true;
                    break;
                }
            }
        }
        if (isAssigned) {
            continue; // Move to the next bond
        }
        // If bond is not assigned, add it to RemainedBonds
        if (!isAssigned) {
            remainedBonds.atomIDs.insert(bond.p1);
            remainedBonds.atomIDs.insert(bond.p2);
            if (!remainedBonds.BondsIndices) {
                remainedBonds.BondsIndices.emplace();
            }
            remainedBonds.BondsIndices->push_back(bondIndex);
        }
    }
}


// Utility function to check if a bond is within a specific range
bool SegmentForceMapper::bondInRange(int p1, int p2, int low, int high) const {
    return (p1 >= low && p1 <= high) && (p2 >= low && p2 <= high);
}

// Utility function to check if a bond is between two specific atoms
bool SegmentForceMapper::bondBetweenAtoms(const BondParams& bond, int atom1, int atom2) const {
    return (bond.p1 == atom1 && bond.p2 == atom2) || (bond.p1 == atom2 && bond.p2 == atom1);
}

// Utility function to allocate bonds to a residue
void SegmentForceMapper::allocateToResidue(PResidues& residue,const int& bondIndex, const BondParams& bond) {
    // Check if bond.p1 or bond.p2 is a hydrogen atom in this residue
    bool p1IsHydrogen = find(residue.HAtomsIDs.begin(), residue.HAtomsIDs.end(), bond.p1) != residue.HAtomsIDs.end();
    bool p2IsHydrogen = find(residue.HAtomsIDs.begin(), residue.HAtomsIDs.end(), bond.p2) != residue.HAtomsIDs.end();


    // Ensure AllBondsIndices is initialized
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
