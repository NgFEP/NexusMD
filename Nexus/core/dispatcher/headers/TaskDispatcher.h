#pragma once

#include <memory>
#include <string>
//#include <vector>
#include <iostream>
#include <SystemXMLParser.h>
#include <PeriodicBoundaryCondition.h>
#include <Coords3D.h>

class TaskDispatcher {
public:
    // Public interface for Simulate with fewer arguments
    virtual void Simulate(const std::string& systemFilename, const std::string& stateFilename,
        const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval) = 0;

    virtual ~TaskDispatcher() = default;

    // Factory method to create TaskDispatcher instances
    static std::shared_ptr<TaskDispatcher> CreateDispatcher(const std::string& processorChoice);

protected:
    // Protected helper methods to initialize default parameters
    virtual void initializeParameters();

    // Simulation parameters stored as members
    std::vector<Coords3D> atomPositions;
    std::vector<double> masses;
    std::vector<PTorsionParams> torsionParams;
    std::vector<BondParams> bondParams;
    std::vector<AngleParams> angleParams;
    NonbondedParams nonbondedParams;
    PeriodicBoundaryCondition::BoxInfo boxInfo;
};
