//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, reporter 

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <limits>
#include "Engine.h"
#include <chrono>


using namespace std;
using namespace BaseLine;;

Engine::Engine(
    const std::string& systemFilename,
    const std::string& stateFilename,
    const std::vector<Coords3D>& atomPositions,
    const std::vector<double>& masses,
    const std::vector<PTorsionParams>& torsionParams,
    const std::vector<BondParams>& bondParams,
    const std::vector<AngleParams>& angleParams,
    const NonbondedParams& nonbondedParams,
    const bool& harmonicBondForceEnabled,
    const bool& harmonicAngleForceEnabled,
    const bool& periodicTorsionForceEnabled,
    const bool& nonbondedForceEnabled
) : _systemFilename(systemFilename),
_stateFilename(stateFilename),
_atomPositions(atomPositions),
_masses(masses),
_torsionParams(torsionParams),
_bondParams(bondParams),
_angleParams(angleParams),
_nonbondedParams(nonbondedParams),
_harmonicBondForceEnabled(harmonicBondForceEnabled),
_harmonicAngleForceEnabled(harmonicAngleForceEnabled),
_periodicTorsionForceEnabled(periodicTorsionForceEnabled),
_nonbondedForceEnabled(nonbondedForceEnabled)
{
    if (_boxInfo.boxSize == Coords3D(0, 0, 0)) { // Check if box size is uninitialized
        // Extract box size from XML if not provided
        //_boxInfo.boxSize = StateXMLParser::extractBoxSize(_stateFilename);
    }
    InitializeSimulationParameters();
}


void Engine::InitializeSimulationParameters() {
    if (!_systemFilename.empty() && !_stateFilename.empty()) {
        _atomPositions = StateXMLParser::parseStateXML(_stateFilename);
        _masses = SystemXMLParser::MassesParser(_systemFilename);
        _torsionParams = SystemXMLParser::PTorsionParser(_systemFilename);
        _bondParams = SystemXMLParser::HBondParser(_systemFilename);// ****
        _angleParams = SystemXMLParser::HAngleParser(_systemFilename);
        _nonbondedParams = SystemXMLParser::NonBondedParser(_systemFilename);
        _RealSimRun = true;

    }
    else {// in this case we're dealing with test files which
        //_boxInfo.boxSize = { 100,100,100 };
    }

    // Extract box boundaries from the initial atom positions
    PeriodicBoundaryCondition::extractBoxBoundaries(_atomPositions, _boxInfo, _stateFilename);

    Initializer initializer;
    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
}

void Engine::CalculateForces() {
    _totalPEnergy = 0.0;


    _totalForces.assign(_numAtoms, Coords3D(0, 0, 0));

    // Forces calculation logic here, updating `totalForces` by adding enabled forces
    if (_harmonicBondForceEnabled) {
        Forces::AddHBond(_totalForces, _atomPositions, _bondParams, _totalPEnergy, _boxInfo);
    }
    if (_harmonicAngleForceEnabled) {
        Forces::AddHAngle(_totalForces, _atomPositions, _angleParams, _totalPEnergy, _boxInfo);
    }
    if (_periodicTorsionForceEnabled) {
        Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams, _totalPEnergy, _boxInfo);
    }
    if (_nonbondedForceEnabled) {
        Forces::AddNonBondElectroPME(_totalForces, _atomPositions, _nonbondedParams, _totalPEnergy, _boxInfo, _exclusions);
    }

}


void Engine::Integrate(int& Step) {
    // Integration logic here, updating `atomPositions` and `velocities`

    VerletIntegration::InverseMasses(_masses, _inverseMasses);
    VerletIntegration::Advance(_atomPositions, _velocities, _totalForces, _inverseMasses, Step, _dt, _boxInfo);

}

void Engine::TotalEnergy(double& timestep) {
    _totalKEnergy = 0.0;
    _totalEnergy = 0.0;
    _kineticEnergies.assign(_numAtoms, 0.0);

    KineticEnergy::calculateKineticEnergy(_velocities, _masses, _totalForces, timestep, _numAtoms, _kineticEnergies, _totalKEnergy);

    _totalEnergy = _totalPEnergy + _totalKEnergy;

}

void Engine::Report(const string& inputFilename, const string& outputFilename, int& step, double& timestep, int& interval) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter reporter;
    if (_RealSimRun) {


        if ((((step + 1) % interval) == 0) || step == 0) {

            //kinetic energy and total energy are calculated only a report of sata is requested at each interval
            TotalEnergy(timestep);
        }

        // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

        if (step == 0) {
            //reporter.pdbOutputGenerator(inputFilename, outputFilename, _outputTemplate, _atomPositions, step);
            reporter.pdbOutputGeneratorPart1(inputFilename, outputFilename, _outputTemplate);
            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
        else if ((((step + 1) % interval) == 0) && step != 0) {
            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }

    }
    else {// if it's a test
        if ((((step + 1) % interval) == 0) || step == 0) {
            //ploting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
    }

}



void Engine::RunSimulation(const string& inputFilename, const string& outputFilename, double& timestep, int& numSteps, int& interval) { // const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs

    // Initialize using the unpacked atomPositions
    _numAtoms = _atomPositions.size();
    _dt = { timestep ,timestep };

    // Initialize _exclusions with the size of _numAtoms
    _exclusions.resize(_numAtoms);

    Exclusions::createExclusions(_numAtoms, _bondParams, _exclusions, _bondCutoff);

    // Loop through each simulation step
    for (int currentStep = 0; currentStep < numSteps; ++currentStep) {

        // Update forces based on current positions
        CalculateForces();

        Report(inputFilename, outputFilename, currentStep, timestep, interval);

        // Update positions and velocities based on new forces
        Integrate(currentStep);

        cout << "Simulation Step #" << currentStep + 1 << " is completed" << endl;


    }
}

