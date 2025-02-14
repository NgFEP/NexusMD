#include "ResidueForceMapper.h"

using namespace std;

ResidueForceMapper::ResidueForceMapper() {};
ResidueForceMapper::~ResidueForceMapper() {};


int ResidueForceMapper::calculateResidueMemory(const Residues& residue) {
    // Size of mandatory vector members
    //int atomMemory = residue.AllAtomsIndices.size();

    //// Total memory for this residue
    //int totalMemory = atomMemory* sizeof(double3);// +angleMemory + torsionMemory;

    //return totalMemory;
    
    
    //// Size of mandatory vector members
    //int atomMemory = residue.AllAtomsIndices.size() * sizeof(int);
    //int hAtomMemory = residue.HAtomsIndices.size() * sizeof(int);
    //int nonHAtomMemory = residue.NonHAtomsIndices.size() * sizeof(int);

    //// Size of optional vector members (check if they have a value before accessing)
    //int bondMemory = residue.AllBondsIndices ? residue.AllBondsIndices->size() * sizeof(int) : 0;
    //int hBondMemory = residue.HBondsIndices ? residue.HBondsIndices->size() * sizeof(int) : 0;
    //int nonHBondMemory = residue.NonHBondsIndices ? residue.NonHBondsIndices->size() * sizeof(int) : 0;
    //int nonResBondAtomMemory = residue.NonResBondAtoms ? residue.NonResBondAtoms->size() * sizeof(int) : 0;

    //int angleMemory = residue.AngleIndices ? residue.AngleIndices->size() * sizeof(int) : 0;
    //int torsionMemory = residue.TorsionIndices ? residue.TorsionIndices->size() * sizeof(int) : 0;

    //// Memory for the string (capacity, since it allocates memory dynamically)
    //int staticMemory = residue.resName.capacity() * sizeof(char);

    //// Total memory for this residue
    //int totalMemory = atomMemory + hAtomMemory + nonHAtomMemory +
    //    bondMemory + hBondMemory + nonHBondMemory +
    //    nonResBondAtomMemory + staticMemory;// +angleMemory + torsionMemory;


    // Size of mandatory vector members
    int atomMemory = residue.AllAtomsIndices.size() * sizeof(double3);
    int hAtomMemory = residue.HAtomsIndices.size() * sizeof(double3);
    int nonHAtomMemory = residue.NonHAtomsIndices.size() * sizeof(double3);

    // Size of optional vector members (check if they have a value before accessing)
    int bondMemory = residue.AllBondsIndices ? residue.AllBondsIndices->size() * sizeof(CudaBondInfo) : 0;
    int bondPEnergyMemory = residue.AllBondsIndices ? residue.AllBondsIndices->size() * sizeof(double) : 0;// in the kernel each bond has its potential energy stored
    int ghgh = sizeof(CudaBondInfo);
    int hBondMemory = residue.HBondsIndices ? residue.HBondsIndices->size() * sizeof(CudaBondInfo) : 0;
    int nonHBondMemory = residue.NonHBondsIndices ? residue.NonHBondsIndices->size() * sizeof(CudaBondInfo) : 0;

    int angleMemory = residue.AngleIndices ? residue.AngleIndices->size() * sizeof(int) : 0;
    int torsionMemory = residue.TorsionIndices ? residue.TorsionIndices->size() * sizeof(int) : 0;


    // Total memory for this residue
    int totalMemory = atomMemory + bondMemory + bondPEnergyMemory;// +angleMemory + torsionMemory;
    int threadSharedMemory = 2 * sizeof(double3) + sizeof(CudaBondInfo) + sizeof(double);

    return totalMemory;
}


void ResidueForceMapper::allocateBonds(const vector<BondParams>& bondParams,vector<Residues>& residues,RemainedBonds& remainedBonds, int& totalResiduesSize, int& totalBondsInResidues, int& totalMemoryOfResidues, vector<vector<int>>& blockResidues) {

    //set<int> assignedConnections;
    set<int> assignedwResidues;
    int numResidues = residues.size();

    //wResidueCounter = 0;
    //connectionCounter = 0;
    //pResidueCounter = 0;
    //bondAtomInResidue = 0;

    for (size_t bondIndex = 0; bondIndex < bondParams.size(); ++bondIndex) {
        const auto& bond = bondParams[bondIndex];
        isAssigned = 0;
        // Check if the bond belongs to any residue
        //if (wResidueCounter==0) {// as soon as water molecules start residue checking ends
        for (size_t i = 0; i < residues.size(); ++i) {
            if (i >= ResAtomIndicesMaxElement.size()) {// each residue has one ResAtomIndicesMaxElement
                ResAtomIndicesMaxElement.push_back(*max_element(residues[i].AllAtomsIndices.begin(), residues[i].AllAtomsIndices.end()));
            }
            if (bond.p1 <= ResAtomIndicesMaxElement[i] || bond.p2 <= ResAtomIndicesMaxElement[i]) {// || ensures connecting bonds are covered in both of the involving residues
                ResBondInfo resBondInfo;
                resBondInfo.d = bond.d;
                resBondInfo.k = bond.k;
                resBondInfo.p1InxGlobal = bond.p1;
                resBondInfo.p2InxGlobal = bond.p2;
                resBondInfo.bondInx = bondIndex;
                bondInResidue(bond.p1, bond.p2, residues[i].AllAtomsIndices, residues[i].NonHAtomsIndices, resBondInfo);
                if (bondAtomInResidue == 2) {
                    allocateToResidue(residues[i], resBondInfo, bond);
                }
                else if (bondAtomInResidue == 1) {
                    if (!residues[i].NonResBondAtoms) {
                        residues[i].NonResBondAtoms.emplace();
                    }
                    residues[i].NonResBondAtoms->push_back(nonResAtomInx);
                    allocateToResidue(residues[i], resBondInfo, bond);
                }
                if (isAssigned == 2) {
                    break;
                }
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
        //if (isAssigned==0 && wResidueCounter < wResidues.size()) {
        //    for (size_t i = wResidueCounter; i < wResidues.size(); ++i) {
        //        if (assignedwResidues.find(i) == assignedwResidues.end() && bondInRangeWResidue(bond.p1, bond.p2, wResidues[i].lowBound, wResidues[i].highBound)) {

        //            if (!wResidues[i].BondsIndices) {
        //                wResidues[i].BondsIndices.emplace();
        //            }

        //            // Store the bond index, not the atom indices
        //            wResidues[i].BondsIndices->push_back(bondIndex);
        //            if (wResidues[i].BondsIndices->size()==2) {//2 bonds alocated to a water molecule
        //                assignedwResidues.insert(i);
        //                wResidueCounter++;
        //            }
        //            isAssigned = 2;
        //            break;
        //        }
        //    }
        //}
        //if (isAssigned==2) {
        //    continue; // Move to the next bond
        //}
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
    sort(residues.begin(), residues.end(), [](const Residues& a, const Residues& b) {
        auto sizeA = a.AllBondsIndices ? a.AllBondsIndices->size() : 0;
        auto sizeB = b.AllBondsIndices ? b.AllBondsIndices->size() : 0;
        return sizeA > sizeB; // Change comparison to greater-than
        });



    for (size_t i = 0; i < residues.size(); ++i) {
        int residueMemory = calculateResidueMemory(residues[i]); // Calculate memory for this residue
        residues[i].resMemSize = residueMemory;
        totalResiduesSize += residueMemory;                // Accumulate total memory
    }


    for (int i = 0; i < numResidues; ++i) {
        totalBondsInResidues += residues[i].AllBondsIndices ? residues[i].AllBondsIndices->size() : 0;
        totalMemoryOfResidues += residues[i].resMemSize ? residues[i].resMemSize : 0;
    }

    int deviceId;
    cudaGetDevice(&deviceId);

    int sharedMemPerBlock;
    cudaDeviceGetAttribute(&sharedMemPerBlock, cudaDevAttrMaxSharedMemoryPerBlock, deviceId);


    // Calculate sizes of static components
    const int sizeOfPositions = sizeof(double3);
    const int sizeOfBondParams = sizeof(BondParams);

    int totalMemoryAllBonds = totalBondsInResidues * (2 * sizeOfPositions + sizeOfBondParams) + totalResiduesSize;
    //int numBlocks = 1.2 * totalMemoryAllBonds / sharedMemPerBlock; // 10% margin
    int numBlocks = 1 + (totalMemoryOfResidues + sharedMemPerBlock - 1) / sharedMemPerBlock;


    // Example of precomputing residue ranges on the host


    //int residueIdx = 0;
    //while(residueIdx < numResidues){
    //    startResidues.push_back(residueIdx);

    //    int sharedMemoryUsed = 0;
    //    while (residueIdx < numResidues && (sharedMemoryUsed + residues[residueIdx].resMemSize < sharedMemPerBlock)) {//0.95 to keep a safe margin of unused sharedMemPerBlock
    //        //sharedMemoryUsed += residues[residueIdx].AllBondsIndices->size() * (2 * sizeOfPositions + sizeOfBondParams) + residues[residueIdx].resMemSize;
    //        sharedMemoryUsed +=  residues[residueIdx].resMemSize;
    //        residueIdx++;
    //    }
    //    endResidues.push_back(residueIdx-1);
    //}

    //vector<vector<int>> blockResidues(numBlocks); // Each block will hold its assigned residues
    blockResidues.resize(numBlocks);

    vector<int> blockSharedMemoryUsed(numBlocks, 0);  // Track shared memory usage for each block

    bool forward = true; // Zigzag direction
    int blockIdx = 0;    // Start with block 0

    for (int residueIdx = 0; residueIdx < numResidues; ++residueIdx) {
        int resMemSize = residues[residueIdx].resMemSize;
        bool allocated = false; // Track if the current residue is allocated

        // Check if adding the residue would exceed the shared memory limit
        if (blockSharedMemoryUsed[blockIdx] + resMemSize < sharedMemPerBlock) {
            blockResidues[blockIdx].push_back(residueIdx);   // Assign residue to the current block
            blockSharedMemoryUsed[blockIdx] += resMemSize;  // Update shared memory usage
            allocated = true; // Mark residue as allocated

        }
        else {
            // If it doesn't fit, try the next block in the current direction
            residueIdx--; // Retry this residue with the next block
        }

        // Update the block index for zigzag pattern
        if (forward) {
            blockIdx++;
            if (blockIdx == numBlocks) { // Reverse direction at the last block
                forward = false;
                blockIdx = numBlocks - 1; // Start moving backward
            }
        }
        else {
            blockIdx--;
            if (blockIdx == -1) { // Reverse direction at the first block
                forward = true;
                blockIdx = 0; // Start moving forward
            }
        }
        // Print message if residue could not be allocated after checking all blocks
        //if (!allocated && blockIdx == 0 && !forward) {
        //    std::cout << "Residue " << residueIdx << " (size: " << resMemSize
        //        << ") could not be allocated to any block. All blocks are full.\n";
        //}
    }

    // Output the block allocations
    //for (int i = 0; i < numBlocks; ++i) {
    //    std::cout << "Block " << i << " residues: ";
    //    for (int resIdx : blockResidues[i]) {
    //        std::cout << resIdx << " ";
    //    }
    //    std::cout << "\nShared memory used: " << blockSharedMemoryUsed[i] << " / " << sharedMemPerBlock << "\n";
    //}







    //cout << ResAtomIndicesMaxElement.size() << "/n" << residues.size();

}


// Utility function to check if a bond is within a specific range of PResidue

//void ResidueForceMapper::bondInResidue(int p1, int p2, const vector<int>& allAtomsIndices) {
//    // Create a hash set for efficient lookup
//    unordered_set<int> atomSet(allAtomsIndices.begin(), allAtomsIndices.end());
//
//    bool p1InResidue = atomSet.find(p1) != atomSet.end();
//    bool p2InResidue = atomSet.find(p2) != atomSet.end();
//
//    if (p1InResidue && p2InResidue) {
//        bondAtomInResidue = 2;
//        isAssigned = 2;
//    }
//    else if (p1InResidue) {
//        bondAtomInResidue = 1;
//        nonResAtomInx = p2;
//        isAssigned++;
//    }
//    else if (p2InResidue) {
//        bondAtomInResidue = 1;
//        nonResAtomInx = p1;
//        isAssigned++;
//    }
//    else {
//        bondAtomInResidue = 0;
//    }
//}




void ResidueForceMapper::bondInResidue(int p1, int p2, vector<int>& allAtomsIndices, vector<int>& nonHAtomsIndices, ResBondInfo& resBondInfo) {
    int p1Index = -1, p2Index = -1; // Initialize indices to -1 (not found)
    bool p1InResidue = false, p2InResidue = false;

    // Single loop to find indices and determine presence
    for (size_t i = 0; i < allAtomsIndices.size(); ++i) {
        if (allAtomsIndices[i] == p1) {
            p1Index = static_cast<int>(i);
            p1InResidue = true;
        }
        if (allAtomsIndices[i] == p2) {
            p2Index = static_cast<int>(i);
            p2InResidue = true;
        }
        // Break early if both are found
        if (p1InResidue && p2InResidue) break;
    }

    // Assign values based on the presence of p1 and p2
    if (p1InResidue && p2InResidue) {
        bondAtomInResidue = 2;
        isAssigned = 2;
        resBondInfo.p1InxLocal = p1Index;
        resBondInfo.p2InxLocal = p2Index;
    }
    else if (p1InResidue) {
        bondAtomInResidue = 1;
        nonResAtomInx = p2;
        isAssigned++;
        resBondInfo.p1InxLocal = p1Index;
        resBondInfo.p2InRes = false;
        resBondInfo.p2InxLocal = allAtomsIndices.size();
        allAtomsIndices.push_back(nonResAtomInx);
        nonHAtomsIndices.push_back(nonResAtomInx); // nonResAtomInx cannot be H


    }
    else if (p2InResidue) {
        bondAtomInResidue = 1;
        nonResAtomInx = p1;
        isAssigned++;
        resBondInfo.p2InxLocal = p2Index;
        resBondInfo.p1InRes = false;
        resBondInfo.p1InxLocal = allAtomsIndices.size();
        allAtomsIndices.push_back(nonResAtomInx);
        nonHAtomsIndices.push_back(nonResAtomInx);

    }
    else {
        bondAtomInResidue = 0;
    }

    // Print or store indices for debugging or further use
    //cout << "p1Index: " << p1Index << ", p2Index: " << p2Index << endl;
}



// Utility function to check if a bond is within a specific range od WResidue
//bool ResidueForceMapper::bondInRangeWResidue(int p1, int p2, int low, int high) {
//    return (p1 >= low && p1 <= high) && (p2 >= low && p2 <= high);
//}

// Utility function to check if a bond is between two specific atoms
//bool ResidueForceMapper::bondBetweenAtoms(const BondParams& bond, int atom1, int atom2) const {
//    return (bond.p1 == atom1 && bond.p2 == atom2) || (bond.p1 == atom2 && bond.p2 == atom1);
//}

// Utility function to allocate bonds to a residue
void ResidueForceMapper::allocateToResidue(Residues& residue,const ResBondInfo& resBondInfo, const BondParams& bond) {

    // Ensure AllBondsIndices is initialized
    if (!residue.AllBondsIndices) {
        //residue.AllBondsIndices = {};
        residue.AllBondsIndices.emplace();

    }
    residue.AllBondsIndices->push_back(resBondInfo);

    // If either atom is a hydrogen atom, classify as a hydrogen bond
    if (residue.resName == "WAT") {
        if (!residue.HBondsIndices) {
            residue.HBondsIndices.emplace();
        }
        residue.HBondsIndices->push_back(resBondInfo);
    }
    else {
        // Check if bond.p1 or bond.p2 is a hydrogen atom in this residue
        bool p1IsHydrogen = find(residue.HAtomsIndices.begin(), residue.HAtomsIndices.end(), bond.p1) != residue.HAtomsIndices.end();
        bool p2IsHydrogen = find(residue.HAtomsIndices.begin(), residue.HAtomsIndices.end(), bond.p2) != residue.HAtomsIndices.end();

        if (p1IsHydrogen || p2IsHydrogen) {
            if (!residue.HBondsIndices) {
                residue.HBondsIndices.emplace();
            }
            residue.HBondsIndices->push_back(resBondInfo);
        }
        else { // Otherwise, classify as a non-hydrogen bond
            if (!residue.NonHBondsIndices) {
                residue.NonHBondsIndices.emplace();
            }
            residue.NonHBondsIndices->push_back(resBondInfo);
        }
    }

}
