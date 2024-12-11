//Engine is the combiner of all the simulation elements: .pdb state reader, initializer(total forces and velocities), total force calculator, integrator, reporter 
// Engine.cpp
//#include "stdafx.h"
#pragma once

#include "CudaEngine.h"

using namespace std;
using namespace Cuda;


//Engine::Engine(const string& SystemFilename, const string& StateFilename) : systemFilename(SystemFilename), stateFilename(StateFilename) {
//    InitializeSimulationParameters();
//}
Engine::Engine(
    const std::string& systemFilename,
    const std::string& stateFilename,
    const std::string& inputFilename,
    const std::vector<Coords3D>& atomPositions,
    const std::vector<double>& masses,
    const std::vector<PTorsionParams>& torsionParams,
    const std::vector<BondParams>& bondParams,
    const std::vector<AngleParams>& angleParams,
    const NonbondedParams& nonbondedParams,
    const PeriodicBoundaryCondition::BoxInfo& boxInfo
) : _systemFilename(systemFilename),
_stateFilename(stateFilename),
_inputFilename(inputFilename),
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
    if (!_systemFilename.empty() && !_stateFilename.empty() && !_inputFilename.empty()) {
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

        // create _atomsBondLoaded
        //_numAtomsBondLoaded = _atomsBondLoaded.size();

        PDBResidueParser pdbResidueParser;
        pdbResidueParser.parseFile(_inputFilename,_residues);// extracting _residues

        _numResidues = _residues.size();
     
        //long startTime = clock();
        auto startTime = chrono::high_resolution_clock::now();
        ResidueForceMapper ResidueForceMapper;
        ResidueForceMapper.allocateBonds(_bondParams, _residues, _remainedBonds, _totalResiduesSize, _totalBondsInResidues, _startResidues, _endResidues);

        //_maxResidueSize = calculateResidueMemory(_residues[0]); // largest residue
        //_maxResidueSize = calculateResidueMemory(_residues.back()); // smallest residue



        //long finishTime = clock();
        auto finishTime = chrono::high_resolution_clock::now();

        // Calculate the elapsed time in microseconds
        auto duration = chrono::duration_cast<chrono::nanoseconds>(finishTime - startTime).count();
        cout << "protein Residue finder Runtime is " << duration << " nanoseconds" << endl;

        //CudaDataStructures CudaDataStructures;
        //CudaDataStructures.bondStructures(_pResidues, _wResidues, _pResiduesBond, _wResiduesBond);

        //// Temporary
        //// verification start
        //// Access the first member of _pResiduesBond
        //const D_PResidues& firstResidue = _pResiduesBond[100];

        //std::cout << "AllBondsIndices of the first member of _pResiduesBond:" << std::endl;

        //// Ensure the pointer is not null and count is valid
        //if (firstResidue.AllBondsIndices != nullptr && firstResidue.AllBondsCount > 0) {
        //    for (int i = 0; i < firstResidue.AllBondsCount; ++i) {
        //        std::cout << "AllBondsIndices[" << i << "] = " << firstResidue.AllBondsIndices[i] << std::endl;
        //    }
        //}
        //// verification end


        std::cout << "Host residue 0 - AllAtomsIndices size: "
            << _residues[0].AllAtomsIndices.size()
            << std::endl;



        CudaBridge CudaBridge;
        cudaMalloc(&d_residues, _numResidues * sizeof(D_Residues));
        CudaBridge.transferResiduesVector(_residues, &d_residues);


        // verification
        //// Define the number of elements to check
        //int numResiduesToCheck = 1; // Adjust as needed

        //// Step 1: Copy the first residue back to the host
        //std::vector<D_Residues> h_residuesCheck(numResiduesToCheck);
        //cudaError_t err = cudaMemcpy(h_residuesCheck.data(), d_residues, numResiduesToCheck * sizeof(D_Residues), cudaMemcpyDeviceToHost);
        //if (err != cudaSuccess) {
        //    std::cerr << "CUDA memcpy failed for residues: " << cudaGetErrorString(err) << std::endl;
        //    cudaFree(d_residues); // Free memory if an error occurs
        //    return;
        //}

        //// Step 2: Validate AllAtomsIndices for the first residue
        //if (h_residuesCheck[0].AllAtomsIndices != nullptr && h_residuesCheck[0].AllAtomsCount > 0) {
        //    // Step 3: Copy AllAtomsIndices data from device to host
        //    std::vector<int> h_AllAtomsIndices(h_residuesCheck[0].AllAtomsCount);
        //    err = cudaMemcpy(h_AllAtomsIndices.data(), h_residuesCheck[0].AllAtomsIndices,
        //        h_residuesCheck[0].AllAtomsCount * sizeof(int), cudaMemcpyDeviceToHost);
        //    if (err != cudaSuccess) {
        //        std::cerr << "CUDA memcpy failed for AllAtomsIndices: " << cudaGetErrorString(err) << std::endl;
        //        return;
        //    }

        //    // Step 4: Output AllAtomsIndices data for verification
        //    std::cout << "AllAtomsIndices for the first residue: ";
        //    for (int i = 0; i < h_residuesCheck[0].AllAtomsCount; ++i) {
        //        std::cout << h_AllAtomsIndices[i] << " ";
        //    }
        //    std::cout << std::endl;

        //    // Step 5: Print the first atom index as a quick check
        //    std::cout << "First atom index in AllAtomsIndices: " << h_AllAtomsIndices[0] << std::endl;
        //}
        //else {
        //    std::cerr << "AllAtomsIndices is nullptr or empty for the first residue." << std::endl;
        //}

        //// Step 6: Print additional details for verification
        //std::cout << "numAtoms in the first residue: " << h_residuesCheck[0].AllAtomsCount << std::endl;
































        // Allocate memory for double3 arrays globally
        _atomPositions_double3 = new double3[_numAtoms];
        _forces_double3 = new double3[_numAtoms];
        _velocities_double3 = new double3[_numAtoms];
        _boxSize_double3 = new double3;
        _lb_double3 = new double3;
        _ub_double3 = new double3;


        // Convert Coords3D (host) to double3 (device)
        for (int i = 0; i < _numAtoms; ++i) {
            _atomPositions_double3[i] = make_double3(_atomPositions[i][0], _atomPositions[i][1], _atomPositions[i][2]);
        }

        // Allocate memory on GPU for the simulation parameters
        cudaMalloc(&d_atomPositions, _numAtoms * sizeof(double3));
        cudaMalloc(&d_masses, _numAtoms * sizeof(double));
        cudaMalloc(&d_inverseMasses, _numAtoms * sizeof(double));
        cudaMalloc(&d_torsionParams, _numTorsions * sizeof(PTorsionParams));  // Assuming TorsionParam is a custom struct
        cudaMalloc(&d_bondParams, _numBonds * sizeof(BondParams));  // Assuming BondParam is a custom struct
        cudaMalloc(&d_angleParams, _numAngles * sizeof(AngleParams));  // Assuming AngleParam is a custom struct
        //cudaMalloc(&d_atomsBondLoaded, _numAtomsBondLoaded * sizeof(ModifiedAtomBondInfo));  // Assuming BondParam is a custom struct
        cudaMalloc(&d_startResidues, _startResidues.size() * sizeof(int));
        cudaMalloc(&d_endResidues, _endResidues.size() * sizeof(int));

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
        cudaMemcpy(d_bondParams, _bondParams.data(), _numBonds * sizeof(BondParams), cudaMemcpyHostToDevice);
        cudaMemcpy(d_angleParams, _angleParams.data(), _numAngles * sizeof(AngleParams), cudaMemcpyHostToDevice);
        //cudaMemcpy(d_atomsBondLoaded, _atomsBondLoaded.data(), _numAtomsBondLoaded * sizeof(ModifiedAtomBondInfo), cudaMemcpyHostToDevice);
        //cudaMemcpy(d_nonbondedParams, _nonbondedParams.data(), numNonbonded * sizeof(NonbondedParam), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numAtoms, &_numAtoms, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numBonds, &_numBonds, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numAngles, &_numAngles, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numTorsions, &_numTorsions, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_numNonbonded, &_numNonbonded, sizeof(int), cudaMemcpyHostToDevice);

        // Copy to device
        cudaMemcpy(d_startResidues, _startResidues.data(), _startResidues.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_endResidues, _endResidues.data(), _endResidues.size() * sizeof(int), cudaMemcpyHostToDevice);



        // instead of copying 
        cudaMemset(d_inverseMasses, 0, _numAtoms * sizeof(double));

        //inverse masses initialization
        launchComputeInverseMassesKernel(d_masses, d_inverseMasses, _numAtoms);

    }
    else {
        // Handle test files case, if any
    }

    // Extract box boundaries from atom positions
    PeriodicBoundaryCondition::extractBoxBoundaries(_atomPositions, _boxInfo, _stateFilename);

    *_boxSize_double3 = make_double3(_boxInfo.boxSize[0], _boxInfo.boxSize[1], _boxInfo.boxSize[2]);
    *_lb_double3 = make_double3(_boxInfo.lb[0], _boxInfo.lb[1], _boxInfo.lb[2]);
    *_ub_double3 = make_double3(_boxInfo.ub[0], _boxInfo.ub[1], _boxInfo.ub[2]);


    // Initialize velocities and forces on CPU (may eventually move to GPU if needed)
    Initializer initializer;
    initializer.InitializeForcesAndVelocities(_atomPositions, _totalForces, _velocities);


    //// Convert Coords3D (forces) to double3 for GPU
    //for (int i = 0; i < _numAtoms; ++i) {
    //    _forces_double3[i] = make_double3(_totalForces[i][0], _totalForces[i][1], _totalForces[i][2]);
    //    _velocities_double3[i] = make_double3(_velocities[i][0], _velocities[i][1], _velocities[i][2]);
    //}

    // Initialize forces and velocities to zero for all atoms
    for (int i = 0; i < _numAtoms; ++i) {
        _forces_double3[i] = make_double3(0.0, 0.0, 0.0);     // Initialize forces to zero
        _velocities_double3[i] = make_double3(0.0, 0.0, 0.0); // Initialize velocities to zero
    }


    // Allocate and copy velocities and forces to GPU
    cudaMalloc(&d_totalForces, _numAtoms * sizeof(double3));
    cudaMalloc(&d_velocities, _numAtoms * sizeof(double3));
    cudaMalloc(&d_kineticEnergies, _numAtoms * sizeof(double));
    cudaMalloc(&d_totalKEnergy, sizeof(double));
    cudaMalloc(&d_totalPEnergy, sizeof(double));
    cudaMalloc(&d_totalEnergy, sizeof(double));



    cudaMalloc(&d_boxsize, sizeof(double3));
    cudaMalloc(&d_lb, sizeof(double3));
    cudaMalloc(&d_ub, sizeof(double3));
    // no need to initialize d_totalForces as it gets initialized every step using cudaMemset in CalculateForces
    //cudaMemcpy(d_totalForces, _forces_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    //cudaMemcpy(d_velocities, _velocities_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_totalPEnergy, &_totalPEnergy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_boxsize, _boxSize_double3, sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lb, _lb_double3, sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ub, _ub_double3, sizeof(double3), cudaMemcpyHostToDevice);


    // instead of copying 
    cudaMemset(d_totalPEnergy, 0, sizeof(double));
    cudaMemset(d_totalKEnergy, 0, sizeof(double));
    cudaMemset(d_totalEnergy, 0, sizeof(double));
    cudaMemset(d_kineticEnergies, 0, _numAtoms * sizeof(double));

    cudaMemset(d_velocities, 0, _numAtoms * sizeof(double3));
    cudaMemset(d_totalForces, 0, _numAtoms * sizeof(double3));

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

    //for (int i = 0; i < _numAtoms; ++i) {
    //    _atomPositions_double3[i] = make_double3(_atomPositions[i][0], _atomPositions[i][1], _atomPositions[i][2]);
    //}
    //// Copy data from host to device
    //cudaMemcpy(d_atomPositions, _atomPositions_double3, _numAtoms * sizeof(double3), cudaMemcpyHostToDevice);


    //long startTime = clock();
    auto startTime = chrono::high_resolution_clock::now();

    // Forces calculation logic here, updating `totalForces` by adding PtorsionForces
    if (!_torsionParams.empty()) {
        //Forces::AddPTorsion(_totalForces, _atomPositions, _torsionParams, _totalPEnergy, _boxInfo);
    }
    if (!_bondParams.empty()) {
        // Forces::AddHBond(_totalForces, _atomPositions, _bondParams, _totalPEnergy, _boxInfo);
        //launchKernelBondForcesGlobal(d_atomPositions, d_bondParams, d_totalForces, d_totalPEnergy,d_boxsize, _numBonds);
        launchKernelBondForcesShared(d_atomPositions, d_bondParams, d_totalForces, d_totalPEnergy, d_boxsize, d_residues, d_startResidues, d_endResidues, _startResidues.size(), _totalBondsInResidues );


    }
    if (!_angleParams.empty()) {
        //Forces::AddHAngle(_totalForces, _atomPositions, _angleParams, _totalPEnergy, _boxInfo);
    }
    if (!_nonbondedParams.particles.empty()) {// if there is no particle infor available to perform the NonBondCalculations
        //Forces::AddNonBondElectroPME(_totalForces, _atomPositions, _nonbondedParams, _totalPEnergy,  _boxInfo,_exclusions);
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
 


    // Function to launch the Verlet integration kernel
    launchVerletIntegrationKernel(d_atomPositions, d_velocities, d_totalForces, d_inverseMasses, d_boxsize, d_lb, d_ub, _numAtoms,  _dt[0]);




    //if (Step == 0) {
    //    VerletIntegration::InverseMasses(_masses, _inverseMasses);
    //}
    //VerletIntegration::Advance(_atomPositions, _velocities, _totalForces, _inverseMasses, Step, _dt, _boxInfo);



}

void Engine::TotalEnergy(double& timestep) {
    _totalKEnergy = 0.0;
    _totalEnergy = 0.0;
    _kineticEnergies.assign(_numAtoms, 0.0);
    cudaMemset(d_totalKEnergy, 0, sizeof(double));


    //KineticEnergy::calculateKineticEnergy(_velocities, _masses, _totalForces, timestep, _numAtoms, _kineticEnergies, _totalKEnergy);

    //_totalEnergy = _totalPEnergy + _totalKEnergy;

    launchKineticEnergyKernel(d_velocities, d_masses, d_inverseMasses, d_totalForces, timestep, _numAtoms, d_kineticEnergies, d_totalKEnergy);

}

void Engine::Report(const string& inputFilename, const string& outputFilename, int& step, double& timestep, int& interval) {
    // Reporting logic here, potentially writing to `outputFilename` for the current `step`
    Reporter reporter;

    //update global parameters for final results reports
    if ((((step + 1) % interval) == 0) || step == 0) {

        //kinetic energy and total energy are calculated only a report of sata is requested at each interval
        TotalEnergy(timestep);

        // Copy results back to host
        cudaMemcpy(_velocities_double3, d_velocities, _numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(_atomPositions_double3, d_atomPositions, _numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(_forces_double3, d_totalForces, _numAtoms * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(&_totalPEnergy, d_totalPEnergy, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&_totalKEnergy, d_totalKEnergy, sizeof(double), cudaMemcpyDeviceToHost);
        //cudaMemcpy(&_kineticEnergies, d_kineticEnergies, _numAtoms * sizeof(double), cudaMemcpyDeviceToHost);

        // Convert double3 results back to Coords3D
        for (size_t i = 0; i < _numAtoms; ++i) {
            _velocities[i] = Coords3D{ _velocities_double3[i].x, _velocities_double3[i].y, _velocities_double3[i].z };
            _atomPositions[i] = Coords3D{ _atomPositions_double3[i].x, _atomPositions_double3[i].y, _atomPositions_double3[i].z };
            _totalForces[i] = Coords3D{ _forces_double3[i].x, _forces_double3[i].y, _forces_double3[i].z };
        }
        _totalEnergy = _totalPEnergy + _totalKEnergy;

    }



    if (_RealSimRun) {

        // Write the REMARK line with the current date
        // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

        if (step == 0) {
            //reporter.pdbOutputGenerator(inputFilename, outputFilename, _outputTemplate, _atomPositions, step);
            reporter.pdbOutputGeneratorPart1(inputFilename, outputFilename, _outputTemplate);
            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
            //plotting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
        else if ((((step + 1) % interval) == 0) && step != 0) {

            _Modelnum = (step + 1) / interval;//in pdb step starts from 1
            reporter.pdbOutputGeneratorPart2(outputFilename, _outputTemplate, _atomPositions, _Modelnum);
            //plotting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
    }
    else {// if it's a test
        if (step == 0) {
            //ploting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
        else if ((((step + 1) % interval) == 0) && step != 0) {
            //ploting energy conservation
            reporter.TotalEnergyReport(outputFilename, _totalKEnergy, _totalPEnergy, _totalEnergy, step);// ******
            //plotting position, velocity and force values
            reporter.TestPVFReport(outputFilename, _atomPositions, _velocities, _totalForces, step, _torsionParams, _bondParams, _angleParams);
        }
    }

}

// Destructor or cleanup function to free GPU memory
void Engine::CleanupGPU() {
    // Free host-side memory if allocated
    if (_atomPositions_double3 != nullptr) {
        delete[] _atomPositions_double3;
        _atomPositions_double3 = nullptr;
    }
    if (_forces_double3 != nullptr) {
        delete[] _forces_double3;
        _forces_double3 = nullptr;
    }
    if (_velocities_double3 != nullptr) {
        delete[] _velocities_double3;
        _velocities_double3 = nullptr;
    }
    if (_boxSize_double3 != nullptr) {
        delete[] _boxSize_double3;
        _boxSize_double3 = nullptr;
    }
    if (_lb_double3 != nullptr) {
        delete[] _lb_double3;
        _lb_double3 = nullptr;
    }
    if (_ub_double3 != nullptr) {
        delete[] _ub_double3;
        _ub_double3 = nullptr;
    }

    // Free GPU memory
    cudaFree(d_atomPositions);
    cudaFree(d_masses);
    cudaFree(d_inverseMasses);
    cudaFree(d_velocities);
    cudaFree(d_totalForces);
    cudaFree(d_bondParams);
    cudaFree(d_torsionParams);
    cudaFree(d_angleParams);
    cudaFree(d_nonbondedParams);
    cudaFree(d_totalPEnergy);
    cudaFree(d_boxsize);
    cudaFree(d_lb);
    cudaFree(d_ub);
    cudaFree(d_numAtoms);
    cudaFree(d_numBonds);
    cudaFree(d_numAngles);
    cudaFree(d_numTorsions);
    cudaFree(d_numNonbonded);
    cudaFree(d_endResidues);
    cudaFree(d_startResidues);

    // Set all GPU pointers to nullptr after freeing
    d_atomPositions = nullptr;
    d_masses = nullptr;
    d_inverseMasses = nullptr;
    d_velocities = nullptr;
    d_totalForces = nullptr;
    d_bondParams = nullptr;
    d_torsionParams = nullptr;
    d_angleParams = nullptr;
    d_nonbondedParams = nullptr;
    d_totalPEnergy = nullptr;
    d_boxsize = nullptr;
    d_lb = nullptr;
    d_ub = nullptr;
    d_numAtoms = nullptr;
    d_numBonds = nullptr;
    d_numAngles = nullptr;
    d_numTorsions = nullptr;
    d_numNonbonded = nullptr;
    d_endResidues = nullptr;
    d_startResidues = nullptr;

    // CUDA Data Structures
    CudaBridge CudaBridge;
    CudaBridge.cleanupPResiduesBond(d_residues, _residues.size());


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

        //long startTime = clock();
        auto startTime2 = chrono::high_resolution_clock::now();

        // Update forces based on current positions
        CalculateForces();
        //TotalEnergy();
        
        // Report current state, clearing the file only at the first step
        // in Report if the interval applies, the _atomPositions, _velocities, _totalPEnergy,_totalKEnergy, _totalEnergy, _totalForces should be updated by CUDA to RAM memory copying 
        Report(inputFilename, outputFilename, currentStep, timestep, interval);

        // Update positions and velocities based on new forces
        //auto [updatedPositions, updatedVelocities] = Integrate(atomPositions, velocities, totalForces, masses, currentStep, timestep);
        Integrate(currentStep);

        // Prepare for the next iteration by updating positions and velocities
        // atomPositions = updatedPositions;
        // velocities = updatedVelocities;

        auto finishTime2 = chrono::high_resolution_clock::now();

        // Calculate the elapsed time in microseconds
        auto duration = chrono::duration_cast<chrono::nanoseconds>(finishTime2 - startTime2).count();
        cout << "1 step simulation Runtime is " << duration << " nanoseconds" << endl;

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
    CleanupGPU();

}

