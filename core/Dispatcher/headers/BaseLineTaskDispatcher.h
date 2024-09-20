#pragma once

#include "TaskDispatcher.h"
#include "Engine.h"
using namespace BaseLine;
class BaseLineTaskDispatcher : public TaskDispatcher {
public:
    void Simulate(const std::string& systemFilename, const std::string& stateFilename, const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval, const std::vector<Coords3D>& atomPositions, const std::vector<double>& masses, const std::vector<PTorsionParams>& torsionParams,
        const std::vector<HBondParams>& bondParams, const std::vector<HAngleParams>& angleParams, const NonbondedParams& nonbondedParams, const PeriodicBoundaryCondition::BoxInfo& boxInfo) override {
        //Engine engine(systemFilename, stateFilename, {}, {}, {}, {}, {});
        Engine engine(systemFilename, stateFilename, atomPositions, masses, torsionParams, bondParams, angleParams, nonbondedParams, boxInfo);


        // Run simulation for CPU
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }
};
