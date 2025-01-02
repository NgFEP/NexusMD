#ifndef CUDADATASTRUCTURES_H
#define CUDADATASTRUCTURES_H

//#include <string>
//#include <vector>
//#include <optional>
//#include "PDBResidueParser.h"
//
//struct ResBondInfo {
//    int bondInx;
//    int p1Inx;
//    int p2Inx;
//};

struct ResBondInfo {
    double d;       // Ideal bond distance in nanometers
    double k;
    int bondInx;
    int p1InxLocal;
    int p2InxLocal;
    int p1InxGlobal;
    int p2InxGlobal;
    bool p1InRes = true;
    bool p2InRes = true;
};

//Unified device Residue Struct for Protein, water, ...
struct D_Residues {
    int* AllAtomsIndices;         // Pointer to all atom indices
    int AllAtomsCount;            // Count of all atoms
    const char* resName;          // 8 bytes (on a 64-bit system)
    ResBondInfo* AllBondsIndices;         // Pointer to all bond indices
    int AllBondsCount;            // Count of all bonds
    int* HAtomsIndices;               // Pointer to Hydrogen atom indices
    int HAtomsCount;              // Count of Hydrogen atoms
    int* NonHAtomsIndices;            // Pointer to Hydrogen atom indices
    int NonHAtomsCount;           // Count of Hydrogen atoms
    ResBondInfo* HBondsIndices;           // Pointer to Hydrogen bond indices
    int HBondsCount;              // Count of Hydrogen bonds
    ResBondInfo* NonHBondsIndices;        // Pointer to non-Hydrogen bond indices
    int NonHBondsCount;           // Count of non-Hydrogen bonds
    int* NonResBondAtoms;         // Pointer to non-residue bond atoms
    int NonResBondAtomsCount;     // Count of non-residue bond atoms
    int resMemSize = 0; // residue memory size

};

struct CudaBondInfo {
    double d;       // Ideal bond distance in nanometers
    double k;
    int bondInx;
    int p1InxLocal;
    int p2InxLocal;
    int p1InxGlobal;
    int p2InxGlobal;
    bool waterMol = false;// does the bond belong to water?
    bool p1InRes = true;
    bool p2InRes = true;
};

//Unified device Residue Struct for Protein, water, ...
struct D_CudaBonds {
    int* AllAtomsIndices;         // Pointer to all atom indices
    int AllAtomsCount;            // Count of all atoms
    CudaBondInfo* AllBondsIndices;         // Pointer to all bond indices
    int AllBondsCount;            // Count of all bonds
    int* HAtomsIndices;               // Pointer to Hydrogen atom indices
    int HAtomsCount;              // Count of Hydrogen atoms
    int* NonHAtomsIndices;            // Pointer to Hydrogen atom indices
    int NonHAtomsCount;           // Count of Hydrogen atoms
    CudaBondInfo* HBondsIndices;           // Pointer to Hydrogen bond indices
    int HBondsCount;              // Count of Hydrogen bonds
    CudaBondInfo* NonHBondsIndices;        // Pointer to non-Hydrogen bond indices
    int NonHBondsCount;           // Count of non-Hydrogen bonds

    // int* NonResBondAtoms;         // Pointer to non-residue bond atoms
    // int NonResBondAtomsCount;     // Count of non-residue bond atoms
    //int resMemSize = 0; // residue memory size

};



//// CUDA-compatible struct for storing residue bond data
//struct D_PResidues {
//    int lowBound;
//    int highBound;
//    const char* resName;          // Pointer to a string (const to indicate it shouldn't be changed on GPU)
//    int* HAtomsIDs;               // Pointer to Hydrogen atom indices
//    int HAtomsCount;              // Count of Hydrogen atoms
//    int* NonHAtomsIDs;               // Pointer to Hydrogen atom indices
//    int NonHAtomsCount;              // Count of Hydrogen atoms
//    int* AllBondsIndices;         // Pointer to all bond indices
//    int AllBondsCount;            // Count of all bonds
//    int* HBondsIndices;           // Pointer to Hydrogen bond indices
//    int HBondsCount;              // Count of Hydrogen bonds
//    int* NonHBondsIndices;        // Pointer to non-Hydrogen bond indices
//    int NonHBondsCount;           // Count of non-Hydrogen bonds
//    int* NonResBondAtoms;         // Pointer to non-residue bond atoms
//    int NonResBondAtomsCount;     // Count of non-residue bond atoms
//};
//
//// CUDA-compatible struct for water residues bond data
//struct D_WResidues {
//    int lowBound;          // Pointer to water residue lower bounds
//    int highBound;         // Pointer to water residue upper bounds
//    //int* HAtomsIDs;               // Pointer to Hydrogen atom indices
//    //int HAtomsCount;              // Count of Hydrogen atoms
//    int* BondsIndices;
//    int BondsCount;
//};

//class CudaDataStructures {
//public:
//    // Function to create CUDA-compatible structures from PResidues and WResidues
//    static void bondStructures(
//        std::vector<PResidues>& pResidues,
//        std::vector<WResidues>& wResidues,
//        std::vector<D_PResidues>& pResiduesBond,
//        std::vector<D_WResidues>& wResiduesBond);
//};

#endif // CUDADATASTRUCTURES_H
