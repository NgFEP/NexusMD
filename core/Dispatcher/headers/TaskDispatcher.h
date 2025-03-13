#pragma once

#include <string>
#include <vector>
#include <memory>
#include <SystemXMLParser.h>
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

    // Functions to enable specific forces
    void EnableHarmonicBondForce() { harmonicBondForceEnabled = true; }
    void EnableHarmonicAngleForce() { harmonicAngleForceEnabled = true; }
    void EnablePeriodicTorsionForce() { periodicTorsionForceEnabled = true; }
    void EnableNonbondedForce() { nonbondedForceEnabled = true; }

    //AndersenThermostat
    void AndersenThermostat(double Temperature = 300, double CollisionFrequency=10) {//kelvin
        _temperature = Temperature;
        _collisionFrequency = CollisionFrequency;
    }

    void MonteCarloBarostat(double Pressure = 1.5, double Frequency = 10) {//bar
        _pressure = Pressure;
        _frequency = Frequency;
    }

protected:
    virtual void initializeParameters();

    // Simulation parameters stored as members
    std::vector<Coords3D> atomPositions;
    std::vector<double> masses;
    std::vector<PTorsionParams> torsionParams;
    std::vector<BondParams> bondParams;
    std::vector<AngleParams> angleParams;
    NonbondedParams nonbondedParams;

    // Individual flags for enabled forces
    bool harmonicBondForceEnabled = false;
    bool harmonicAngleForceEnabled = false;
    bool periodicTorsionForceEnabled = false;
    bool nonbondedForceEnabled = false;
    double _temperature;
    double _collisionFrequency;
    double _pressure;
    double _frequency;

};
