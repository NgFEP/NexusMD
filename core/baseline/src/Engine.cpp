//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, reporter 
// Engine.cpp
#include "stdafx.h"
#include "Engine.h"

using namespace std;
using namespace BaseLine;;


//Engine::Engine(const string& SystemFilename, const string& StateFilename) : systemFilename(SystemFilename), stateFilename(StateFilename) {
//    InitializeSimulationParameters();
//}
Engine::Engine(
    const string& systemFilename,
    const string& stateFilename,
    const vector<Coords3D>& atomPositions,
    const vector<double>& masses,
    const vector<PTorsionParams>& torsionParams,
    const vector<HAngleParams>& angleParams
) : _systemFilename(systemFilename),
_stateFilename(stateFilename),
_atomPositions(atomPositions),
_masses(masses),
_torsionParams(torsionParams),
_angleParams(angleParams)
{
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
        _angleParams = SystemXMLParser::HAngleParser(_systemFilename);
    }
    Initializer initializer;
    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
}



//void Engine::Initialize(vector<Coords3D>& atomPositions, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {
//    // Directly call Initializer::InitializeForcesAndVelocities without needing a return value.
//    // This directly updates totalForces and velocities.
//    Initializer initializer;
//
//    initializer.InitializeForcesAndVelocities(atomPositions, totalForces, velocities);
//}


void Engine::CalculateForces() {
    
    //bug found: Forces are calculated from scratch every step and we need to make sure that the vector gets equal to zero on each step before adding up new calculated forces.
    size_t numAtoms = _atomPositions.size();
    _totalForces.assign(numAtoms, Coords3D(0, 0, 0));
    
    // Forces calculation logic here, updating `totalForces` by adding PtorsionForces
    if (!_torsionParams.empty()) {
        Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams);
    }
    if (!_angleParams.empty()) {
        Forces::AddHAngle(_totalForces, _atomPositions, _angleParams);
    }
}

//pair<vector<Coords3D>, vector<Coords3D>> Engine::Integrate2(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& masses, int& Step, double& StepSize) {
//    // Integration logic here, updating `atomPositions` and `velocities`
//    VerletIntegration::advance(atomPositions, velocities, totalForces, masses, Step, StepSize);
//    return pair(atomPositions, velocities);
//}

void Engine::Integrate(int& Step, double& StepSize) {
    // Integration logic here, updating `atomPositions` and `velocities`
    VerletIntegration::advance(_atomPositions, _velocities, _totalForces, _masses, Step, StepSize);
}

void Engine::Report(const string& outputFilename, int step) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter reporter;
    reporter.report(outputFilename, _atomPositions, _velocities, _totalForces,step, _torsionParams, _angleParams);

}



void Engine::RunSimulation(const string& outputFilename, double timestep, int numSteps){ //, const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs

    // Initialize using the unpacked atomPositions
    //Initialize(atomPositions, totalForces, velocities);

    // Loop through each simulation step
    for (int currentStep = 0; currentStep < numSteps; ++currentStep) {

        // Update forces based on current positions
        CalculateForces();

        // Report current state, clearing the file only at the first step
        Report(outputFilename, currentStep);

        // Update positions and velocities based on new forces
        //auto [updatedPositions, updatedVelocities] = Integrate(atomPositions, velocities, totalForces, masses, currentStep, timestep);
        Integrate(currentStep, timestep);


        // Prepare for the next iteration by updating positions and velocities
        //atomPositions = updatedPositions;
        //velocities = updatedVelocities;

    }
}

