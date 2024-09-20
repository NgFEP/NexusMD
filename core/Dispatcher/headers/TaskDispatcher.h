#pragma once

#include <memory>
#include <string>
#include <iostream>
#include <Engine.h>
//using namespace BaseLine;
class TaskDispatcher {
public:
    virtual void Simulate(const std::string& systemFilename, const std::string& stateFilename,
        const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval,const std::vector<Coords3D>& atomPositions, const std::vector<double>& masses, const std::vector<PTorsionParams>& torsionParams,
        const std::vector<HBondParams>& bondParams, const std::vector<HAngleParams>& angleParams , const NonbondedParams& nonbondedParams, const PeriodicBoundaryCondition::BoxInfo& boxInfo) = 0;

    virtual ~TaskDispatcher() = default;

    static std::unique_ptr<TaskDispatcher> CreateDispatcher(const std::string& processorChoice);
};
