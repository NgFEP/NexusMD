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
#include "Coords3D.h"
#include "Reporter.h"

namespace BaseLine {
    class Engine {
    public:
        // Constructor declaration
        Engine(const std::string& systemFilename, const std::string& stateFilename);
        tuple<std::vector<Coords3D>, std::vector<double>, std::vector<TorsionParameters>> Inputs(const std::string& systemFilename, const std::string& stateFilename);

        //void Inputs(const std::string& systemFilename, const std::string& stateFilename);
        //void RunSimulation(const std::string& outputFilename, double timestep, int numSteps, const std::string& systemFilename, const std::string& stateFilename);
        void RunSimulation(const string& outputFilename, double timestep, int numSteps, const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities);


    private:
        std::string systemFilename;
        std::string stateFilename;
        std::vector<Coords3D> atomPositions;
        std::vector<TorsionParameters> torsionParams;
        std::vector<Coords3D> totalForces;
        std::vector<Coords3D> velocities;
        std::vector<double> masses;
        //void Initialize(const std::vector<Coords3D>& atomPositions);
        //std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Initialize(const std::vector<Coords3D>& atomPositions);
        void Initialize(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);

        std::vector<Coords3D> CalculateForces(const vector<Coords3D>& atomPositions, vector<TorsionParameters>& torsionParams, vector<Coords3D>& totalForces);
        //std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Integrate(int& StepNum, double& StepSize);
        std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Integrate(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, vector<double>& masses, int& Step, double& StepSize);
        //void Report(const std::string& outputFilename, int step);
        void Report(const std::string& outputFilename, std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, int step, vector<TorsionParameters>& torsionParams);

    };
}

#endif // ENGINE_H
