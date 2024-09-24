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
#include <chrono>
#include "CudaEngine.h"

using namespace std;
using namespace Cuda;;


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
        //_boxInfo.boxSize = StateXMLParser::extractBoxSize(_stateFilename);
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




//void Engine::InitializeSimulationParameters() {
//    if (!_systemFilename.empty() && !_stateFilename.empty()) {
//        _atomPositions = StateXMLParser::parseStateXML(_stateFilename);
//        _masses = SystemXMLParser::MassesParser(_systemFilename);
//        _torsionParams = SystemXMLParser::PTorsionParser(_systemFilename);
//        _bondParams = SystemXMLParser::HBondParser(_systemFilename);// ****
//        _angleParams = SystemXMLParser::HAngleParser(_systemFilename);
//        _nonbondedParams = SystemXMLParser::NonBondedParser(_systemFilename);
//        _RealSimRun = true;
//
//    }
//    else {// in this case we're dealing with test files which
//        //_boxInfo.boxSize = { 100,100,100 };
//    }
//
//    // Extract box boundaries from the initial atom positions
//    PeriodicBoundaryCondition::extractBoxBoundaries(_atomPositions, _boxInfo, _stateFilename);
//
//    Initializer initializer;
//    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);
//}


void Engine::InitializeSimulationParameters() {
    if (!_systemFilename.empty() && !_stateFilename.empty()) {
        // Parse data from input files
        _atomPositions = StateXMLParser::parseStateXML(_stateFilename);
        _masses = SystemXMLParser::MassesParser(_systemFilename);
        _torsionParams = SystemXMLParser::PTorsionParser(_systemFilename);
        _bondParams = SystemXMLParser::HBondParser(_systemFilename);
        _angleParams = SystemXMLParser::HAngleParser(_systemFilename);
        _nonbondedParams = SystemXMLParser::NonBondedParser(_systemFilename);
        _RealSimRun = true;

        _numAtoms = _atomPositions.size();
        _numBonds = _bondParams.size();
        _numAngles = _angleParams.size();
        _numTorsions = _torsionParams.size();
        _numNonbonded = _nonbondedParams.particles.size();

        // Allocate memory for double3 arrays globally
        _atomPositions_double3 = new double3[_numAtoms];
        _forces_double3 = new double3[_numAtoms];
        _velocities_double3 = new double3[_numAtoms];
        _boxSize_double3 = new double3;

        // Convert Coords3D (host) to double3 (device)
        for (int i = 0; i < _numAtoms; ++i) {
            _atomPositions_double3[i] = make_double3(_atomPositions[i][0], _atomPositions[i][1], _atomPositions[i][2]);
        }

        // Allocate memory on GPU for the simulation parameters
        cudaMalloc(&d_atomPositions, _numAtoms * sizeof(double3));
        cudaMalloc(&d_masses, _numAtoms * sizeof(double));
        cudaMalloc(&d_torsionParams, _numTorsions * sizeof(PTorsionParams));  // Assuming TorsionParam is a custom struct
        cudaMalloc(&d_bondParams, _numBonds * sizeof(HBondParams));  // Assuming BondParam is a custom struct
        cudaMalloc(&d_angleParams, _numAngles * sizeof(HAngleParams));  // Assuming AngleParam is a custom struct
        //cudaMalloc(&d_nonbondedParams, numNonbonded * sizeof(NonbondedParam));  // Assuming NonbondedParam is a custom struct
        cudaMalloc(&d_numAtoms, sizeof(int));
        cudaMalloc(&d_numBonds, sizeof(int));
        cudaMalloc(&d_numAngles, sizeof(int));
        cudaMalloc(&d_numTorsions, sizeof(int));
        cudaMalloc(&d_numNonbonded, sizeof(int));


        // Copy the simulation parameters from CPU to GPU
        cudaMemcpy(d_atomPositions, _atomPositions_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
        cudaMemcpy(d_masses, _masses.data(), _numAtoms * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_torsionParams, _torsionParams.data(), _numTorsions * sizeof(PTorsionParams), cudaMemcpyHostToDevice);
        cudaMemcpy(d_bondParams, _bondParams.data(), _numBonds * sizeof(HBondParams), cudaMemcpyHostToDevice);
        cudaMemcpy(d_angleParams, _angleParams.data(), _numAngles * sizeof(HAngleParams), cudaMemcpyHostToDevice);
        //cudaMemcpy(d_nonbondedParams, _nonbondedParams.data(), numNonbonded * sizeof(NonbondedParam), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numAtoms, &_numAtoms, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numBonds, &_numBonds, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numAngles, &_numAngles, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numTorsions, &_numTorsions, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numNonbonded, &_numNonbonded, sizeof(int), cudaMemcpyHostToDevice);


    }
    else {
        // Handle test files case, if any
    }

    // Extract box boundaries from atom positions
    PeriodicBoundaryCondition::extractBoxBoundaries(_atomPositions, _boxInfo, _stateFilename);

    *_boxSize_double3 = make_double3(_boxInfo.boxSize[0], _boxInfo.boxSize[1], _boxInfo.boxSize[2]);


    // Initialize velocities and forces on CPU (may eventually move to GPU if needed)
    Initializer initializer;
    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);


    // Convert Coords3D (forces) to double3 for GPU
    for (int i = 0; i < _numAtoms; ++i) {
        _forces_double3[i] = make_double3(_totalForces[i][0], _totalForces[i][1], _totalForces[i][2]);
        _velocities_double3[i] = make_double3(_velocities[i][0], _velocities[i][1], _velocities[i][2]);
    }

    // Allocate and copy velocities and forces to GPU
    cudaMalloc(&d_totalForces, _numAtoms * sizeof(double3));
    cudaMalloc(&d_velocities, _numAtoms * sizeof(double3));
    cudaMalloc(&d_totalPEnergy, sizeof(double));
    cudaMalloc(&d_boxsize, sizeof(double3));
    // no need to initialize d_totalForces as it gets initialized every step using cudaMemset in CalculateForces
    //cudaMemcpy(d_totalForces, _forces_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocities, _velocities_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_totalPEnergy, &_totalPEnergy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_boxsize, _boxSize_double3, sizeof(double3), cudaMemcpyHostToDevice);



    // Free allocated memory for global arrays on the host
    //delete[] _atomPositions_double3;
    //delete[] _forces_double3;
    //delete[] _velocities_double3;
    //delete _boxSize_double3;

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
    cudaMemset(d_totalPEnergy, 0, sizeof(double));
    //bug found: Forces are calculated from scratch every step and we need to make sure that the vector gets equal to zero on each step before adding up new calculated forces.
    _totalForces.assign(_numAtoms, Coords3D(0, 0, 0));
    cudaMemset(d_totalForces, 0, _numAtoms * sizeof(double3));


    //everystep d_atomPositions should be updated since the data comes from verletintegration and it's still on CPU

    for (int i = 0; i < _numAtoms; ++i) {
        _atomPositions_double3[i] = make_double3(_atomPositions[i][0], _atomPositions[i][1], _atomPositions[i][2]);
    }
    // Copy data from host to device
    cudaMemcpy(d_atomPositions, _atomPositions_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);




    //long startTime = clock();
    auto startTime = chrono::high_resolution_clock::now();

    // Forces calculation logic here, updating `totalForces` by adding PtorsionForces
    if (!_torsionParams.empty()) {
        //Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams, _totalPEnergy, _boxInfo);
    }
    if (!_bondParams.empty()) {
        // Forces::AddHBond(_totalForces, _atomPositions, _bondParams, _totalPEnergy, _boxInfo);

        launchKernelBondForces_V2(d_atomPositions, d_bondParams, d_totalForces, d_totalPEnergy,d_boxsize, _numBonds);

    }
    if (!_angleParams.empty()) {
        //Forces::AddHAngle(_totalForces, _atomPositions, _angleParams, _totalPEnergy, _boxInfo);
    }
    if (!_nonbondedParams.particles.empty()) {// if there is no particle infor available to perform the NonBondCalculations
        //Forces::AddNonBondElectroPME(_totalForces, _atomPositions, _nonbondedParams, _totalPEnergy,  _boxInfo,_exclusions);
    }

    //for now _totalForces and _totalPEnergy should be updated every step since kinetic energy and verlet integration are in cpu and need the updated values
    // if the interval applies, the _atomPositions, _velocities, _totalPEnergy,_totalKEnergy, _totalEnergy, _totalForces should be updated 
    // Copy results back to host
    cudaMemcpy(_forces_double3, d_totalForces, _numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
    cudaMemcpy(&_totalPEnergy, d_totalPEnergy, sizeof(double), cudaMemcpyDeviceToHost);

    // Convert double3 results back to Coords3D
    for (size_t i = 0; i < _numAtoms; ++i) {
        _totalForces[i] = Coords3D{ _forces_double3[i].x, _forces_double3[i].y, _forces_double3[i].z };
    }



    //long finishTime = clock();
    auto finishTime = chrono::high_resolution_clock::now();

    // Calculate the elapsed time in microseconds
    auto duration = chrono::duration_cast<chrono::nanoseconds>(finishTime - startTime).count();
    cout << "Force Calculation Runtime is " << duration << " nanoseconds" << endl;

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
    if (_RealSimRun) {

        // Write the REMARK line with the current date
        // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

        if (step == 0) {

            //reporter.pdbOutputGenerator(inputFilename, outputFilename, _outputTemplate, _atomPositions, step);
            reporter.pdbOutputGeneratorPart1(inputFilename, outputFilename, _outputTemplate);
            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);

            //ploting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******

            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
        else if ((((step + 1) % interval) == 0) && step != 0) {

            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
            //ploting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******

            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
    }
}



void Engine::RunSimulation(const string& inputFilename, const string& outputFilename, double& timestep, int& numSteps, int& interval){ // const string& systemFilename, const string& stateFilename, vector<Coords3D>& totalForces, vector<Coords3D>& velocities) {    // Unpack the tuple returned by Inputs

    // Initialize using the unpacked atomPositions
    // Initialize(atomPositions, totalForces, velocities);
    //_numAtoms = _atomPositions.size();
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
        // in Report if the interval applies, the _atomPositions, _velocities, _totalPEnergy,_totalKEnergy, _totalEnergy, _totalForces should be updated by CUDA to RAM memory copying 
        Report(inputFilename, outputFilename, currentStep, interval);

        // Update positions and velocities based on new forces
        //auto [updatedPositions, updatedVelocities] = Integrate(atomPositions, velocities, totalForces, masses, currentStep, timestep);
        Integrate(currentStep);


        // Prepare for the next iteration by updating positions and velocities
        // atomPositions = updatedPositions;
        // velocities = updatedVelocities;

    }

    // Synchronize the device to wait for all kernel operations to complete
    cudaDeviceSynchronize();

    // Check for any errors returned by the kernel launch
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        // Print the error message if any CUDA error occurred
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
    }

    // Free device memory
    cudaFree(d_atomPositions);
    cudaFree(d_bondParams);
    cudaFree(d_totalForces);
    cudaFree(d_totalPEnergy);
    cudaFree(d_boxsize);

}

