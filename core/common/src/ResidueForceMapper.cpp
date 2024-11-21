#include "ResidueForceMapper.h"

using namespace std;

ResidueForceMapper::ResidueForceMapper() {};
ResidueForceMapper::~ResidueForceMapper() {};


void ResidueForceMapper::allocateBonds(const vector<BondParams>& bondParams,vector<PResidues>& pResidues,vector<WResidues>& wResidues,RemainedBonds& remainedBonds) {

    set<int> assignedConnections;
    set<int> assignedwResidues;

    wResidueCounter = 0;
    //connectionCounter = 0;
    //pResidueCounter = 0;
    //bondAtomInRangeCounter = 0;


    for (size_t bondIndex = 0; bondIndex < bondParams.size(); ++bondIndex) {
        const auto& bond = bondParams[bondIndex];
        isAssigned = 0;

        // Check if the bond belongs to any residue
        //if (wResidueCounter==0) {// as soon as water molecules start residue checking ends
        for (size_t i = 0; i < pResidues.size(); ++i) {
            bondInRangePResidue(bond.p1, bond.p2, pResidues[i].lowBound, pResidues[i].highBound);
            if (bondAtomInRangeCounter == 2) {
                allocateToResidue(pResidues[i], bondIndex, bond);
            }
            else if (bondAtomInRangeCounter == 1) {
                if (!pResidues[i].NonResBondAtoms) {
                    pResidues[i].NonResBondAtoms.emplace();
                }
                pResidues[i].NonResBondAtoms->push_back(nonResAtomInx);
                allocateToResidue(pResidues[i], bondIndex, bond);

            }
            if (isAssigned == 2) {
                break;
            }
        }
        //}
        if (isAssigned==2) {
            continue; // Move to the next bond
        }
        //// Check if the bond belongs to any connection
        //if (isAssigned==0 && connectionCounter < connections.size()) {
        //    for (size_t i = connectionCounter; i < connections.size(); ++i) {
        //        if (assignedConnections.find(i) == assignedConnections.end() &&
        //            bondBetweenAtoms(bond, connections[i].atomC, connections[i].atomN)) {
        //            // Store the bond index, not the atom indices
        //            connections[i].BondIndex = bondIndex;
        //            assignedConnections.insert(i);
        //            connectionCounter++;
        //            //pResidueCounter++;
        //            isAssigned = true;
        //            break;
        //        }
        //    }
        //}
        //if (isAssigned==2) {
        //    continue; // Move to the next bond
        //}
        // Check if the bond belongs to any water molecule
        if (isAssigned==0 && wResidueCounter < wResidues.size()) {
            for (size_t i = wResidueCounter; i < wResidues.size(); ++i) {
                if (assignedwResidues.find(i) == assignedwResidues.end() && bondInRangeWResidue(bond.p1, bond.p2, wResidues[i].lowBound, wResidues[i].highBound)) {

                    if (!wResidues[i].BondsIndices) {
                        wResidues[i].BondsIndices.emplace();
                    }

                    // Store the bond index, not the atom indices
                    wResidues[i].BondsIndices->push_back(bondIndex);
                    if (wResidues[i].BondsIndices->size()==2) {//2 bonds alocated to a water molecule
                        assignedwResidues.insert(i);
                        wResidueCounter++;
                    }
                    isAssigned = 2;
                    break;
                }
            }
        }
        if (isAssigned==2) {
            continue; // Move to the next bond
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
}


// Utility function to check if a bond is within a specific range of PResidue
void ResidueForceMapper::bondInRangePResidue(int p1, int p2, int low, int high) {
    if ((p1 >= low && p1 <= high) && (p2 >= low && p2 <= high)) {
        bondAtomInRangeCounter = 2;
        isAssigned = 2;
    }
    else if ((p1 >= low && p1 <= high)) {
        bondAtomInRangeCounter = 1;
        nonResAtomInx = p2;
        isAssigned++;
    }
    else if ((p2 >= low && p2 <= high)) {
        bondAtomInRangeCounter = 1;
        nonResAtomInx = p1;
        isAssigned++;
    }
    else {
        bondAtomInRangeCounter = 0;
    }
}

// Utility function to check if a bond is within a specific range od WResidue
bool ResidueForceMapper::bondInRangeWResidue(int p1, int p2, int low, int high) {
    return (p1 >= low && p1 <= high) && (p2 >= low && p2 <= high);
}

// Utility function to check if a bond is between two specific atoms
bool ResidueForceMapper::bondBetweenAtoms(const BondParams& bond, int atom1, int atom2) const {
    return (bond.p1 == atom1 && bond.p2 == atom2) || (bond.p1 == atom2 && bond.p2 == atom1);
}

// Utility function to allocate bonds to a residue
void ResidueForceMapper::allocateToResidue(PResidues& pResidue,const int& bondIndex, const BondParams& bond) {
    // Check if bond.p1 or bond.p2 is a hydrogen atom in this residue
    bool p1IsHydrogen = find(pResidue.HAtomsIDs.begin(), pResidue.HAtomsIDs.end(), bond.p1) != pResidue.HAtomsIDs.end();
    bool p2IsHydrogen = find(pResidue.HAtomsIDs.begin(), pResidue.HAtomsIDs.end(), bond.p2) != pResidue.HAtomsIDs.end();


    // Ensure AllBondsIndices is initialized
    if (!pResidue.AllBondsIndices) {
        //residue.AllBondsIndices = {};
        pResidue.AllBondsIndices.emplace();

    }
    pResidue.AllBondsIndices->push_back(bondIndex);

    // If either atom is a hydrogen atom, classify as a hydrogen bond
    if (p1IsHydrogen || p2IsHydrogen) {
        if (!pResidue.HBondsIndices) {
            pResidue.HBondsIndices.emplace();
        }
        pResidue.HBondsIndices->push_back(bondIndex);
    }

    else { // Otherwise, classify as a non-hydrogen bond
        if (!pResidue.NonHBondsIndices) {
            pResidue.NonHBondsIndices.emplace();
        }
        pResidue.NonHBondsIndices->push_back(bondIndex);
    }
}
