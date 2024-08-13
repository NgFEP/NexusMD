//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, reporter 
// Engine.cpp
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

using namespace std;
using namespace BaseLine;;


//Engine::Engine(const string& SystemFilename, const string& StateFilename) : systemFilename(SystemFilename), stateFilename(StateFilename) {
//    InitializeSimulationParameters();
//}
Engine::Engine(
    const std::string& systemFilename,
    const std::string& stateFilename,
    const std::vector<Coords3D>& atomPositions,
    const std::vector<double>& masses,
    const std::vector<PTorsionParams>& torsionParams,
    const std::vector<HBondParams>& bondParams,
    const std::vector<HAngleParams>& angleParams,
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
        _boxInfo.boxSize = StateXMLParser::extractBoxSize(_stateFilename);
    }
    InitializeSimulationParameters();
}


//void Engine::InitializeFromFiles(const string& systemFilename, const string& stateFilename) {
//    _atomPositions = StateXMLParser::parseStateXML(stateFilename);
//    _masses = SystemXMLParser::MassesParser(systemFilename);
//    _torsionParams = SystemXMLParser::PTorsionParser(systemFilename);
//    _angleParams = SystemXMLParser::HAngleParser(systemFilename);
//    Initializer initializer;
//    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
//}



//tuple<vector<Coords3D>, vector<double>, vector<PTorsionParams>> Engine::Inputs(const string& systemFilename, const string& stateFilename) {
//    // Assuming StateXMLParser and SystemXMLParser are classes that have been correctly implemented and their methods are static
//    auto atomPositions = StateXMLParser::parseStateXML(stateFilename);
//    auto masses = SystemXMLParser::MassesParser(systemFilename);
//    auto torsionParams = SystemXMLParser::PTorsionParser(systemFilename);
//
//    return make_tuple(atomPositions, masses, torsionParams);
//}




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



//void Engine::InitializeSimulationParameters() {
//    if (!_systemFilename.empty() && !_stateFilename.empty()) {
//        _atomPositions = StateXMLParser::parseStateXML(_stateFilename);
//        _masses = SystemXMLParser::MassesParser(_systemFilename);
//        _torsionParams = SystemXMLParser::PTorsionParser(_systemFilename);
//        _angleParams = SystemXMLParser::HAngleParser(_systemFilename);
//        _RealSimRun = true;
//       
//    }
//    Initializer initializer;
//    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
//}



//void Engine::Initialize(vector<Coords3D>& atomPositions, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {
//    // Directly call Initializer::InitializeForcesAndVelocities without needing a return value.
//    // This directly updates totalForces and velocities.
//    Initializer initializer;
//
//    initializer.InitializeForcesAndVelocities(atomPositions, totalForces, velocities);
//}


void Engine::CalculateForces() {
    //since energy is calculated during force calculation for each step _totalPEnergy gets reseted
    //potential energy is calculated per force and per atom energy is not calculated therefore only totalPEnergy is going to be reported
    _totalPEnergy = 0.0;


    //bug found: Forces are calculated from scratch every step and we need to make sure that the vector gets equal to zero on each step before adding up new calculated forces.
    _totalForces.assign(_numAtoms, Coords3D(0, 0, 0));
    
    // Forces calculation logic here, updating `totalForces` by adding PtorsionForces
    if (!_torsionParams.empty()) {
        // Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams, _totalPEnergy, _boxInfo);
    }
    if (!_bondParams.empty()) {
        // Forces::AddHBond(_totalForces, _atomPositions, _bondParams, _totalPEnergy, _boxInfo);
    }
    if (!_angleParams.empty()) {
        // Forces::AddHAngle(_totalForces, _atomPositions, _angleParams, _totalPEnergy, _boxInfo);
    }
    if (!_nonbondedParams.particles.empty()) {// if there is no particle infor available to perform the NonBondCalculations
        Forces::AddNonBondElectroPME(_totalForces, _atomPositions, _nonbondedParams, _totalPEnergy,  _boxInfo,_exclusions);
    }



}

//pair<vector<Coords3D>, vector<Coords3D>> Engine::Integrate2(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& masses, int& Step, double& StepSize) {
//    // Integration logic here, updating `atomPositions` and `velocities`
//    VerletIntegration::advance(atomPositions, velocities, totalForces, masses, Step, StepSize);
//    return pair(atomPositions, velocities);
//}

void Engine::Integrate(int& Step) {
    // Integration logic here, updating `atomPositions` and `velocities`

    VerletIntegration::InverseMasses(_masses, _inverseMasses);
    //VerletIntegration::Advance(_dt, _atomPositions, _velocities, _totalForces, _inverseMasses, _posDelta);
    VerletIntegration::Advance(_atomPositions, _velocities, _totalForces, _inverseMasses, Step, _dt, _boxInfo);
    //VerletIntegration::SelectVerletStepSize(_velocities, _totalForces, _inverseMasses, _dt,  errorTol,  maxStepSize);

   // VerletIntegration::advance(_atomPositions, _velocities, _totalForces, _masses, Step, StepSize);
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
    //vector<double> _kineticEnergies;// a vector of _kineticEnergy for each atom
    //vector<double> _potentialEnergies;// a vector of _potentialEnergy for each atom

}

void Engine::Report(const string& inputFilename, const string& outputFilename, int& step, int& interval) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter reporter;
    if (_RealSimRun){

        // Write the REMARK line with the current date
        // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

        if (step == 0 ) {
            //reporter.pdbOutputGenerator(inputFilename, outputFilename, _outputTemplate, _atomPositions, step);

            reporter.pdbOutputGeneratorPart1(inputFilename, outputFilename, _outputTemplate);
        }
        else if (((step + 1) % interval) == 0) {
            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
        }
    }
    else {
        reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
    }
    //ploting energy conservation
    reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******

}



void Engine::RunSimulation(const string& inputFilename, const string& outputFilename, double& timestep, int& numSteps, int& interval){ // const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs

    // Initialize using the unpacked atomPositions
    // Initialize(atomPositions, totalForces, velocities);
    _numAtoms = _atomPositions.size();
    _dt = { timestep ,timestep };

    // Initialize _exclusions with the size of _numAtoms
    _exclusions.resize(_numAtoms);

    // Exclusions::createExclusions(_numAtoms, _bondParams, _angleParams, _torsionParams, _exclusions, _bondCutoff);
    Exclusions::createExclusions(_numAtoms, _bondParams, _exclusions, _bondCutoff);

    // Loop through each simulation step
    for (int currentStep = 0; currentStep < numSteps; ++currentStep) {

        // Update forces based on current positions
        CalculateForces();
        //TotalEnergy();
        TotalEnergy(timestep);

        // Report current state, clearing the file only at the first step
        
        Report(inputFilename, outputFilename, currentStep, interval);

        // Update positions and velocities based on new forces
        //auto [updatedPositions, updatedVelocities] = Integrate(atomPositions, velocities, totalForces, masses, currentStep, timestep);
        Integrate(currentStep);


        // Prepare for the next iteration by updating positions and velocities
        // atomPositions = updatedPositions;
        // velocities = updatedVelocities;

    }
}

