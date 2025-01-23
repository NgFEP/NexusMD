#ifndef CUDA_BRIDGE_H
#define CUDA_BRIDGE_H

#include <vector>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cuda_runtime.h>
#include "PDBResidueParser.h"
#include "CudaDataStructures.h"
#include "CudaForceMapper.h"

namespace Cuda {

    class CudaBridge {
    public:
        // Constructor and Destructor
        CudaBridge();
        ~CudaBridge();

        // Functions
        void transferResiduesVector(const std::vector<Residues>& hostResidues, D_Residues** deviceResiduesArray);
        void transferCudaBondsVector(const std::vector<CudaBonds>& h_bonds, D_CudaBonds** d_bondsArray);
        void cleanupPResiduesBond(D_Residues* deviceResiduesArray, size_t numResidues);

    private:
        template <typename T>
        void copyVectorToDevice(const std::vector<T>& hostVector, T** devicePointer);
        void transferSingleResidue(const Residues& hostResidue, D_Residues& deviceResidue);
        void transferSingleCudaBond(const CudaBonds& h_bond, D_CudaBonds& d_bond);

    };
}

#endif // CUDA_BRIDGE_H
