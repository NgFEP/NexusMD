//#pragma once
//
//#include "TaskDispatcher.h"
//#include "Engine.h"
//
//class CudaTaskDispatcher : public TaskDispatcher {
//public:
//    void Simulate(const std::string& systemFilename, const std::string& stateFilename, const std::string& inputFilename, const std::string& outputFilename,
//        double StepSize, int TotalSteps, int interval) override {
//        BaseLine::Engine engine(systemFilename, stateFilename, {}, {}, {}, {}, {});
//
//        // Run simulation for CPU
//        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
//    }
//};
//
//






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

        // Use the initialized parameters in Cuda::Engine
        Cuda::Engine engine(systemFilename, stateFilename, inputFilename, atomPositions,
            masses, torsionParams, bondParams, angleParams,
            nonbondedParams, boxInfo);

        // Run simulation for CUDA
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }

protected:
    // Customize parameter initialization for CudaTaskDispatcher
    void initializeParameters() override {
        // Initialize parameters specific to CUDA
        atomPositions = {};  // Initialize as needed
        masses = {};         // Initialize as needed
        torsionParams = {};  // Initialize as needed
        bondParams = {};     // Initialize as needed
        angleParams = {};    // Initialize as needed
        nonbondedParams = {}; // Initialize as needed
        boxInfo.boxSize = {}; // Initialize as needed
    }
};
