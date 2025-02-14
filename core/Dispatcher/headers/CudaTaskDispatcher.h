#pragma once

#include "TaskDispatcher.h"
#include "CudaEngine.h"

class CudaTaskDispatcher : public TaskDispatcher {
public:
    // Override the simplified Simulate method
    void Simulate(const std::string& systemFilename, const std::string& stateFilename,
        const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval) override {
        // Initialize parameters
        initializeParameters();


        if (harmonicBondForceEnabled) {
            std::cout << "Harmonic Bond Force is enabled" << std::endl;
        }
        if (harmonicAngleForceEnabled) {
            std::cout << "Harmonic Angle Force is enabled" << std::endl;
        }
        if (periodicTorsionForceEnabled) {
            std::cout << "Periodic Torsion Force is enabled" << std::endl;
        }
        if (nonbondedForceEnabled) {
            std::cout << "Nonbonded Force is enabled" << std::endl;
        }


        // Use the initialized parameters in Cuda::Engine
        Cuda::Engine engine(systemFilename, stateFilename, inputFilename, atomPositions,
            masses, torsionParams, bondParams, angleParams,
            nonbondedParams, harmonicBondForceEnabled,
            harmonicAngleForceEnabled, periodicTorsionForceEnabled, nonbondedForceEnabled);


        // Run simulation for CUDA
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }

protected:
    // Customize parameter initialization for CudaTaskDispatcher
    void initializeParameters() override {

        atomPositions = {};  // Initialize as needed
        masses = {};         // Initialize as needed
        torsionParams = {};  // Initialize as needed
        bondParams = {};     // Initialize as needed
        angleParams = {};    // Initialize as needed
        nonbondedParams = {}; // Initialize as needed
    }
};
