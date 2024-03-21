//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, reporter 
// Engine.cpp
#include "Engine.h"

using namespace std;
using namespace BaseLine;;


Engine::Engine(const string& SystemFilename, const string& StateFilename) : systemFilename(SystemFilename), stateFilename(StateFilename) {}



tuple<vector<Coords3D>, vector<double>, vector<TorsionParameters>> Engine::Inputs(const string& systemFilename, const string& stateFilename) {
    // Assuming StateXMLParser and SystemXMLParser are classes that have been correctly implemented and their methods are static
    auto atomPositions = StateXMLParser::parseStateXML(stateFilename);
    auto masses = SystemXMLParser::MassesParser(systemFilename);
    auto torsionParams = SystemXMLParser::PeriodicTorsionParser(systemFilename);

    return make_tuple(atomPositions, masses, torsionParams);
}




void Engine::Initialize(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities) {
    // Directly call Initializer::InitializeForcesAndVelocities without needing a return value.
    // This directly updates totalForces and velocities.
    Initializer initializer;

    initializer.InitializeForcesAndVelocities(atomPositions, totalForces, velocities);
}


vector<Coords3D> Engine::CalculateForces(const vector<Coords3D>& atomPositions, vector<TorsionParameters>& torsionParams, vector<Coords3D>& totalForces) {
    // Forces calculation logic here, updating `totalForces`
    Forces::AddPeriodicTorsion(totalForces, atomPositions, torsionParams);

    return totalForces;
}

pair<vector<Coords3D>, vector<Coords3D>> Engine::Integrate(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& masses, int& Step, double& StepSize) {
    // Integration logic here, updating `atomPositions` and `velocities`
    VerletIntegration::advance(atomPositions, velocities, totalForces, masses, Step, StepSize);
    return pair(atomPositions, velocities);
}

void Engine::Report(const string& outputFilename, vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, int step, bool clearFile) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter::report(outputFilename, atomPositions, velocities, totalForces,step, clearFile);

}



void Engine::RunSimulation(const string& outputFilename, double timestep, int numSteps, const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs
    auto [atomPositions, masses, torsionParams] = Inputs(systemFilename, stateFilename);

    // Initialize using the unpacked atomPositions
    Initialize(atomPositions, totalForces, velocities);

    // Loop through each simulation step
    for (int currentStep = 0; currentStep < numSteps; ++currentStep) {
        // Update forces based on current positions
        CalculateForces(atomPositions, torsionParams, totalForces);

        // Update positions and velocities based on new forces
        auto [updatedPositions, updatedVelocities] = Integrate(atomPositions, velocities, totalForces, masses, currentStep, timestep);

        // Report current state, clearing the file only at the first step
        bool clearFile = (currentStep == 0);
        Report(outputFilename, updatedPositions, updatedVelocities, totalForces, currentStep, clearFile);

        // Prepare for the next iteration by updating positions and velocities
        atomPositions = updatedPositions;
        velocities = updatedVelocities;
    }
}

