#ifndef CUDA_FORCE_MAPPER_H
#define CUDA_FORCE_MAPPER_H

#include "PDBResidueParser.h" // Include for residues
#include <vector>


// Unified host CudaBonds Struct for storing bond information of residues inside the blocks
struct CudaBonds {
    std::vector<int> AllAtomsIndices; // List of atom indices
    std::vector<int> HAtomsIndices; // List of Hydrogen atom indices
    std::vector<int> NonHAtomsIndices; // List of Non Hydrogen atom indices
    std::vector<CudaBondInfo> AllBondsIndices; // List of all bond indices
    //std::optional<std::vector<int>> AllBondsIndices; // List of all bond indices
    std::vector<CudaBondInfo> HBondsIndices; // List of Hydrogen bond indices
    std::vector<CudaBondInfo> NonHBondsIndices;// List of non-Hydrogen bond indices
};

class CudaForceMapper {
public:
    // Constructor
    CudaForceMapper();
    ~CudaForceMapper();

    void cudaAllocateBonds(std::vector<Residues>& residues, std::vector<CudaBonds>& cudaBonds, std::vector<int>& startResidues, std::vector<int>& endResidues);


private:

};

#endif // CUDA_FORCE_MAPPER_H
