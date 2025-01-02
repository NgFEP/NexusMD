#include "CudaBridge.h"

using namespace std;
using namespace Cuda;

// Constructor
CudaBridge::CudaBridge() {
    // Initialization if needed
}

// Destructor
CudaBridge::~CudaBridge() {
    // Cleanup if needed
}

// Helper function to copy vector data from host to device
template <typename T>
void CudaBridge::copyVectorToDevice(const vector<T>& hostVector, T** devicePointer) {
    // Check if the input vector is not empty
    if (!hostVector.empty()) {
        // Allocate memory on the device
        cudaError_t err = cudaMalloc(reinterpret_cast<void**>(devicePointer), hostVector.size() * sizeof(T));
        if (err != cudaSuccess) {
            cerr << "CUDA malloc failed: " << cudaGetErrorString(err) << endl;
            *devicePointer = nullptr; // Ensure the device pointer is null to avoid undefined behavior
            return;
        }

        // Copy data from host to device
        err = cudaMemcpy(*devicePointer, hostVector.data(), hostVector.size() * sizeof(T), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            cerr << "CUDA memcpy failed: " << cudaGetErrorString(err) << endl;
            cudaFree(*devicePointer); // Free allocated memory on failure
            *devicePointer = nullptr; // Set pointer to null
        }
    }
    else {
        // If the input vector is empty, set the device pointer to nullptr
        *devicePointer = nullptr;
    }
}


// Transfer a single Residue to D_Residues
void CudaBridge::transferSingleResidue(const Residues& h_residue, D_Residues& d_residue) {

    d_residue.AllAtomsIndices = nullptr;
    copyVectorToDevice(h_residue.AllAtomsIndices, &d_residue.AllAtomsIndices);
    d_residue.AllAtomsCount = h_residue.AllAtomsIndices.size();

    // Transfer string
    d_residue.resName = nullptr;
    cudaMalloc((void**)&d_residue.resName, h_residue.resName.size() + 1);
    cudaMemcpy((void*)d_residue.resName, h_residue.resName.c_str(), h_residue.resName.size() + 1, cudaMemcpyHostToDevice);

    // Transfer vectors
    d_residue.HAtomsIndices = nullptr;
    copyVectorToDevice(h_residue.HAtomsIndices, &d_residue.HAtomsIndices);
    d_residue.HAtomsCount = h_residue.HAtomsIndices.size();

    d_residue.NonHAtomsIndices = nullptr;
    copyVectorToDevice(h_residue.NonHAtomsIndices, &d_residue.NonHAtomsIndices);
    d_residue.NonHAtomsCount = h_residue.NonHAtomsIndices.size();

    if (h_residue.AllBondsIndices.has_value()) {
        copyVectorToDevice(*h_residue.AllBondsIndices, &d_residue.AllBondsIndices);
        d_residue.AllBondsCount = h_residue.AllBondsIndices->size();
    }
    else {
        d_residue.AllBondsIndices = nullptr;
        d_residue.AllBondsCount = 0;
    }
    if (h_residue.HBondsIndices.has_value()) {
        copyVectorToDevice(*h_residue.HBondsIndices, &d_residue.HBondsIndices);
        d_residue.HBondsCount = h_residue.HBondsIndices->size();
    }
    else {
        d_residue.HBondsIndices = nullptr;
        d_residue.HBondsCount = 0;
    }

    if (h_residue.NonHBondsIndices.has_value()) {
        copyVectorToDevice(*h_residue.NonHBondsIndices, &d_residue.NonHBondsIndices);
        d_residue.NonHBondsCount = h_residue.NonHBondsIndices->size();
    }
    else {
        d_residue.NonHBondsIndices = nullptr;
        d_residue.NonHBondsCount = 0;
    }

    if (h_residue.NonResBondAtoms.has_value()) {
        copyVectorToDevice(*h_residue.NonResBondAtoms, &d_residue.NonResBondAtoms);
        d_residue.NonResBondAtomsCount = h_residue.NonResBondAtoms->size();
    }
    else {
        d_residue.NonResBondAtoms = nullptr;
        d_residue.NonResBondAtomsCount = 0;
    }

    //d_residue.resMemSize = nullptr;
    d_residue.resMemSize = h_residue.resMemSize;

}

void CudaBridge::transferResiduesVector(const vector<Residues>& h_residues, D_Residues** d_residuesArray) {
    int numResidues = h_residues.size();
    vector<D_Residues> d_residues(numResidues);

    for (size_t i = 0; i < numResidues; ++i) {
        transferSingleResidue(h_residues[i], d_residues[i]);
    }

    // Allocate device array for D_Residues
    cudaError_t err = cudaMalloc(d_residuesArray, numResidues * sizeof(D_Residues));
    if (err != cudaSuccess) {
        cerr << "CUDA malloc failed for D_Residues: " << cudaGetErrorString(err) << endl;
        return;
    }

    // Copy the array of D_Residues to device
    err = cudaMemcpy(*d_residuesArray, d_residues.data(), numResidues * sizeof(D_Residues), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cerr << "CUDA memcpy failed: " << cudaGetErrorString(err) << endl;
    }
}

// Transfer a single CudaBonds to D_CudaBonds
void CudaBridge::transferSingleCudaBond(const CudaBonds& h_bond, D_CudaBonds& d_bond) {
    // Transfer vectors
    d_bond.AllAtomsIndices = nullptr;
    copyVectorToDevice(h_bond.AllAtomsIndices, &d_bond.AllAtomsIndices);
    d_bond.AllAtomsCount = h_bond.AllAtomsIndices.size();

    d_bond.HAtomsIndices = nullptr;
    copyVectorToDevice(h_bond.HAtomsIndices, &d_bond.HAtomsIndices);
    d_bond.HAtomsCount = h_bond.HAtomsIndices.size();

    d_bond.NonHAtomsIndices = nullptr;
    copyVectorToDevice(h_bond.NonHAtomsIndices, &d_bond.NonHAtomsIndices);
    d_bond.NonHAtomsCount = h_bond.NonHAtomsIndices.size();

    // Transfer bond information arrays
    d_bond.AllBondsIndices = nullptr;
    copyVectorToDevice(h_bond.AllBondsIndices, &d_bond.AllBondsIndices);
    d_bond.AllBondsCount = h_bond.AllBondsIndices.size();

    d_bond.HBondsIndices = nullptr;
    copyVectorToDevice(h_bond.HBondsIndices, &d_bond.HBondsIndices);
    d_bond.HBondsCount = h_bond.HBondsIndices.size();

    d_bond.NonHBondsIndices = nullptr;
    copyVectorToDevice(h_bond.NonHBondsIndices, &d_bond.NonHBondsIndices);
    d_bond.NonHBondsCount = h_bond.NonHBondsIndices.size();
}

// Transfer the vector of CudaBonds to the device
void CudaBridge::transferCudaBondsVector(const vector<CudaBonds>& h_cudaBonds, D_CudaBonds** d_cudaBondsArray) {
    int numCudaBonds = h_cudaBonds.size();
    vector<D_CudaBonds> d_bonds(numCudaBonds);

    // Transfer each bond individually
    for (size_t i = 0; i < numCudaBonds; ++i) {
        transferSingleCudaBond(h_cudaBonds[i], d_bonds[i]);
    }

    // Allocate memory for the device array of D_CudaBonds
    cudaError_t err = cudaMalloc(d_cudaBondsArray, numCudaBonds * sizeof(D_CudaBonds));
    if (err != cudaSuccess) {
        cerr << "CUDA malloc failed for D_CudaBonds: " << cudaGetErrorString(err) << endl;
        return;
    }

    // Copy the array of D_CudaBonds to the device
    err = cudaMemcpy(*d_cudaBondsArray, d_bonds.data(), numCudaBonds * sizeof(D_CudaBonds), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        cerr << "CUDA memcpy failed: " << cudaGetErrorString(err) << endl;
    }
}




// Cleanup function
void CudaBridge::cleanupPResiduesBond(D_Residues* d_residuesArray, size_t numResidues) {
    //vector<Residues> d_pResiduesBond(numResidues);

    //// Copy the device array to host
    //cudaMemcpy(d_pResiduesBond.data(), d_residuesArray, numResidues * sizeof(Residues), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < numResidues; ++i) {
        if (d_residuesArray[i].resName) cudaFree((void*)d_residuesArray[i].resName);
        if (d_residuesArray[i].AllAtomsIndices) cudaFree(d_residuesArray[i].AllAtomsIndices);
        if (d_residuesArray[i].HAtomsIndices) cudaFree(d_residuesArray[i].HAtomsIndices);
        if (d_residuesArray[i].NonHAtomsIndices) cudaFree(d_residuesArray[i].NonHAtomsIndices);
        if (d_residuesArray[i].AllBondsIndices) cudaFree(d_residuesArray[i].AllBondsIndices);
        if (d_residuesArray[i].HBondsIndices) cudaFree(d_residuesArray[i].HBondsIndices);
        if (d_residuesArray[i].NonHBondsIndices) cudaFree(d_residuesArray[i].NonHBondsIndices);
        if (d_residuesArray[i].NonResBondAtoms) cudaFree(d_residuesArray[i].NonResBondAtoms);
        //if (d_residuesArray[i].resMemSize) cudaFree((void*)d_residuesArray[i].resMemSize);

    }

    cudaFree(d_residuesArray);
}

