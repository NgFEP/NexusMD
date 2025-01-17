
#pragma once

#include "TaskDispatcher.h"
#include "Engine.h"

using namespace BaseLine;

class BaseLineTaskDispatcher : public TaskDispatcher {
public:
    // To override the simplified Simulate method
    void Simulate(const std::string& systemFilename, const std::string& stateFilename,
        const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval) override {
        // To initialize parameters
        initializeParameters();

        // Use the initialized parameters in Engine
        Engine engine(systemFilename, stateFilename, atomPositions, masses, torsionParams,
            bondParams, angleParams, nonbondedParams, boxInfo);

        // To run simulation for CPU
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }

protected:
    // To customize parameter initialization for BaseLineTaskDispatcher
    void initializeParameters() override {
        // To initialize parameters specific to BaseLine
        atomPositions = {};  // Initialize as needed
        masses = {};         // Initialize as needed
        torsionParams = {};  // Initialize as needed
        bondParams = {};     // Initialize as needed
        angleParams = {};    // Initialize as needed
        nonbondedParams = {}; // Initialize as needed
        boxInfo.boxSize = {}; // Initialize as needed
    }
};
