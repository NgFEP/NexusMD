//#pragma once
//
//#include "TaskDispatcher.h"
//#include "Engine.h"
//
//class CudaTaskDispatcher : public TaskDispatcher {
//public:
//    void Simulate(const std::string& systemFilename, const std::string& stateFilename, const std::string& inputFilename, const std::string& outputFilename,
//        double StepSize, int TotalSteps, int interval) override {
//        BaseLine::Engine engine(systemFilename, stateFilename, {}, {}, {}, {}, {});
//
//        // Run simulation for CPU
//        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
//    }
//};
//
//






#pragma once

#include "TaskDispatcher.h"
#include "CudaEngine.h"
#include "Coords3D.h"
//using namespace Cuda;
class CudaTaskDispatcher : public TaskDispatcher {
public:
    void Simulate(const std::string& systemFilename, const std::string& stateFilename, const std::string& inputFilename, const std::string& outputFilename,
        double StepSize, int TotalSteps, int interval, const std::vector<Coords3D>& atomPositions, const std::vector<double>& masses, const std::vector<PTorsionParams>& torsionParams, 
        const std::vector<HBondParams>& bondParams, const std::vector<HAngleParams>& angleParams , const NonbondedParams& nonbondedParams, const PeriodicBoundaryCondition::BoxInfo& boxInfo) override {
        Cuda::Engine engine(systemFilename, stateFilename, atomPositions, masses, torsionParams, bondParams, angleParams, nonbondedParams, boxInfo);
        
        // Run simulation for CUDA
        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);
    }
};
