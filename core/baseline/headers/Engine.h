//Engine is the heart of the NexaBind which does the actual simulation (evolving the system by combining the OpenMM system data, accumulating forces, and verlet integrator)
// Engine.h
#ifndef ENGINE_H
#define ENGINE_H

#include "StateXMLParser.h"
#include "SystemXMLParser.h"
#include "PeriodicTorsionForce.h"
#include "Initializer.h"
#include "Forces.h"
#include "VerletIntegration.h"
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include "Coords3D.h"
#include "Reporter.h"

namespace BaseLine {
    class Engine {
    public:
        // Constructor declaration
        //Engine(const std::string& systemFilename, const std::string& stateFilename);
        Engine(//this decleration allows optional inputs for engine object construction useful for both test and production runs
            const std::string& systemFilename = "",
            const std::string& stateFilename = "",
            const std::vector<Coords3D>& atomPositions = std::vector<Coords3D>(),
            const std::vector<double>& masses = std::vector<double>(),
            const std::vector<PTorsionParams>& torsionParams = std::vector<PTorsionParams>(),
            const std::vector<HAngleParams>& angleParams = std::vector<HAngleParams>()
        );

        
        
        //void Inputs(const string& systemFilename, const string& stateFilename);
        //void RunSimulation(const string& outputFilename, double timestep, int numSteps, const string& systemFilename, const string& stateFilename);
        void RunSimulation(const std::string& outputFilename, double timestep, int numSteps);// , const std::string& systemFilename, const std::string& stateFilename, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);


    private:
        void InitializeSimulationParameters();
        std::string _systemFilename;
        std::string _stateFilename;
        std::vector<Coords3D> _atomPositions;
        std::vector<PTorsionParams> _torsionParams;
        std::vector<HAngleParams> _angleParams;
        std::vector<Coords3D> _totalForces;
        std::vector<Coords3D> _velocities;
        std::vector<double> _masses;
        //void InitializeFromFiles(const std::string& systemFilename, const std::string& stateFilename);

        //void Initialize(const vector<Coords3D>& atomPositions);
        //pair<vector<Coords3D>, vector<Coords3D>> Initialize(const vector<Coords3D>& atomPositions);
        //void Initialize(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
        //std::vector<Coords3D> CalculateForces(const std::vector<Coords3D>& atomPositions, const std::vector<PTorsionParams>& torsionParams, const std::vector<HAngleParams>& angleParams, std::vector<Coords3D>& totalForces);
        void CalculateForces();
        //pair<vector<Coords3D>, vector<Coords3D>> Integrate(int& StepNum, double& StepSize);
        //std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Integrate2(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& Step, double& StepSize);
        void Engine::Integrate(int& Step, double& StepSize);
        //void Report(const string& outputFilename, int step);
        void Report(const std::string& outputFilename, int step);

    };
}

#endif // ENGINE_H
