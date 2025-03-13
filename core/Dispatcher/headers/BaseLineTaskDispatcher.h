
#pragma once

#include "TaskDispatcher.h"
#include "Engine.h"

using namespace BaseLine;

class BaseLineTaskDispatcher : public TaskDispatcher {
public:
    // Override
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

        // Using initialized parameters in Engine
        Engine engine(systemFilename, stateFilename, inputFilename, atomPositions, masses, torsionParams,
            bondParams, angleParams, nonbondedParams, harmonicBondForceEnabled,
            harmonicAngleForceEnabled, periodicTorsionForceEnabled, nonbondedForceEnabled, _temperature, _collisionFrequency, _pressure, _frequency);

        // Run simulation for CPU
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }

protected:
    // To customize parameter initialization for BaseLineTaskDispatcher
    void initializeParameters() override {
        atomPositions = {}; 
        masses = {};         
        torsionParams = {};
        bondParams = {};  
        angleParams = {};   
        nonbondedParams = {}; 
    }
};
