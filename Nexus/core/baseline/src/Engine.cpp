//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, and reporter 
//#include "stdafx.h"
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
    const PeriodicBoundaryCondition::BoxInfo& boxInfo
) : _systemFilename(systemFilename),
_stateFilename(stateFilename),
_atomPositions(atomPositions),
_masses(masses),
_torsionParams(torsionParams),
_bondParams(bondParams),
_angleParams(angleParams),
_nonbondedParams(nonbondedParams),
_boxInfo(boxInfo)
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

    // To extract box boundaries from the initial atom positions
    PeriodicBoundaryCondition::extractBoxBoundaries(_atomPositions, _boxInfo, _stateFilename);

    Initializer initializer;
    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
}



void Engine::CalculateForces() {
    //since energy is calculated during force calculation for each step _totalPEnergy gets reseted
    //potential energy is calculated per force and per atom energy is not calculated therefore only totalPEnergy is going to be reported
    _totalPEnergy = 0.0;


    //To reset Forces
    _totalForces.assign(_numAtoms, Coords3D(0, 0, 0));

    //long startTime = clock();
    auto startTime = chrono::high_resolution_clock::now();

    // Forces calculation logic here, updating `totalForces` by adding PtorsionForces
    if (!_torsionParams.empty()) {
        //Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams, _totalPEnergy, _boxInfo);
    }
    if (!_bondParams.empty()) {
        Forces::AddHBond(_totalForces, _atomPositions, _bondParams, _totalPEnergy, _boxInfo);
    }
    if (!_angleParams.empty()) {
        //Forces::AddHAngle(_totalForces, _atomPositions, _angleParams, _totalPEnergy, _boxInfo);
    }
    if (!_nonbondedParams.particles.empty()) {// if there is no particle infor available to perform the NonBondCalculations
        //Forces::AddNonBondElectroPME(_totalForces, _atomPositions, _nonbondedParams, _totalPEnergy, _boxInfo, _exclusions);
    }

    //long finishTime = clock();
    auto finishTime = chrono::high_resolution_clock::now();

    // To calculate the elapsed time in microseconds
    auto duration = chrono::duration_cast<chrono::nanoseconds>(finishTime - startTime).count();
    cout << "Force Calculation Runtime is " << duration << " nanoseconds" << endl;


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

    //KineticEnergy::calculateKineticEnergy(_velocities, _masses,_numAtoms, _kineticEnergies, _totalKEnergy);
    KineticEnergy::calculateKineticEnergy(_velocities, _masses, _totalForces, timestep, _numAtoms, _kineticEnergies, _totalKEnergy);

    _totalEnergy = _totalPEnergy + _totalKEnergy;
    if (_totalEnergy > 200) {
        cout << "";
    }

}

void Engine::Report(const string& inputFilename, const string& outputFilename, int& step,double& timestep, int& interval) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter reporter;
    if (_RealSimRun){


        if ((((step + 1) % interval) == 0) || step == 0) {

            //kinetic energy and total energy are calculated only a report of sata is requested at each interval
            TotalEnergy(timestep);
        }

        // To write the REMARK line with the current date
        // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

        if (step == 0 ) {
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



void Engine::RunSimulation(const string& inputFilename, const string& outputFilename, double& timestep, int& numSteps, int& interval){ // const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs

    // To initialize using the unpacked atomPositions
    _numAtoms = _atomPositions.size();
    _dt = { timestep ,timestep };

    // To initialize _exclusions with the size of _numAtoms
    _exclusions.resize(_numAtoms);

    Exclusions::createExclusions(_numAtoms, _bondParams, _exclusions, _bondCutoff);

    // To loop through each simulation step
    for (int currentStep = 0; currentStep < numSteps; ++currentStep) {


        //long startTime = clock();
        auto startTime2 = chrono::high_resolution_clock::now();


        // To update forces based on current positions
        CalculateForces();
        //TotalEnergy();

        // To report current state, clearing the file only at the first step
        
        Report(inputFilename, outputFilename, currentStep, timestep, interval);

        // To update positions and velocities based on new forces
        Integrate(currentStep);


        // To prepare for the next iteration by updating positions and velocities

        auto finishTime2 = chrono::high_resolution_clock::now();
        // To calculate the elapsed time in microseconds
        auto duration = chrono::duration_cast<chrono::nanoseconds>(finishTime2 - startTime2).count();
        cout << "1 step simulation Runtime is " << duration << " nanoseconds" << endl;

    }
}

