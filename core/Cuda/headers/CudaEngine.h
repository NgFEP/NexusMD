//Engine is the heart of the NexaBind which does the actual simulation (evolving the system by combining the OpenMM system data, accumulating forces, and verlet integrator)
// Engine.h
#ifndef CUDAENGINE_H
#define CUDAENGINE_H

#include "CudaStateXMLParser.h"
#include "SystemXMLParser.h"
#include "PDBResidueParser.h"
#include "ResidueForceMapper.h"
#include "CudaForceMapper.h"
#include "CudaDataStructures.h"
#include "CudaBridge.h"
#include "CudaPeriodicTorsionForce.h"
#include "CudaInitializer.h"
#include "CudaForces.h"
#include "PeriodicBoundaryCondition.h"
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <chrono>
#include "Coords3D.h"
#include "CudaReporter.h"
#include "CudaExclusions.h"
#include "InverseMassesKernel.h"
//#include "CudaModifiedForceDataStructures.h"
#include "HarmonicBondForceKernel.h"
#include "VerletIntegrationKernel.h"
#include "KineticEnergyKernel.h"


namespace Cuda {

    //struct BoxInfo {
    //    Coords3D lb;// lower boundary
    //    Coords3D ub;// upper boundary
    //    Coords3D boxSize;

    //    BoxInfo() : lb(0.0, 0.0, 0.0), ub(0.0, 0.0, 0.0), boxSize(0.0, 0.0, 0.0) {}
    //    BoxInfo(const Coords3D& lower, const Coords3D& upper)
    //        : lb(lower), ub(upper), boxSize(std::abs(upper[0] - lower[0]), std::abs(upper[1] - lower[1]), std::abs(upper[2] - lower[2])) {}
    //};



    class Engine {
    public:
        // Constructor declaration
        //Engine(const std::string& systemFilename, const std::string& stateFilename);
        Engine(//this decleration allows optional inputs for engine object construction useful for both test and production runs
            const std::string& systemFilename = "",
            const std::string& stateFilename = "",
            const std::string& inputFile = "",
            const std::vector<Coords3D>& atomPositions = std::vector<Coords3D>(),
            const std::vector<double>& masses = std::vector<double>(),
            const std::vector<PTorsionParams>& torsionParams = std::vector<PTorsionParams>(),
            const std::vector<BondParams>& bondParams = std::vector<BondParams>(),
            const std::vector<AngleParams>& angleParams = std::vector<AngleParams>(),
            const NonbondedParams& nonbondedParams = NonbondedParams(),
            bool harmonicBondForceEnabled = false,
            bool harmonicAngleForceEnabled = false,
            bool periodicTorsionForceEnabled = false,
            bool nonbondedForceEnabled = false

        );
        
        
        //void Inputs(const string& systemFilename, const string& stateFilename);
        //void RunSimulation(const string& outputFilename, double timestep, int numSteps, const string& systemFilename, const string& stateFilename);
        void RunSimulation(const std::string& inputFilename, const std::string& outputFilename, double& timestep, int& numSteps, int& interval);// , const std::string& systemFilename, const std::string& stateFilename, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);


    private:
        // CPU data members

        void InitializeSimulationParameters();
        bool _RealSimRun=false;
        int _numAtoms;
        int _numAtomsBondLoaded;
        int _numBonds;
        int _numAngles;
        int _numTorsions;
        int _numNonbonded;
        int _numResidues;
        //int _maxResidueSize;
        int _totalBondsInResidues = 0;
        int _totalMemoryOfResidues = 0;

        int _totalResiduesSize = 0;

        std::string _systemFilename;
        std::string _stateFilename;
        std::string _inputFilename;
        std::vector<Coords3D> _atomPositions;
        std::vector<PTorsionParams> _torsionParams;
        std::vector<BondParams> _bondParams;
        std::vector<AngleParams> _angleParams;
        NonbondedParams _nonbondedParams;
        std::vector<std::set<int>> _exclusions;
        //std::vector<double> _bondPEnergies;

        // PDBResidueParser
        std::vector<Residues> _residues; 
        std::vector<CudaBonds> _cudaBonds;
        // ResidueForceMapper
        RemainedBonds _remainedBonds;
        // CudaDatastructures
        //std::vector<int> _startResidues;
        //std::vector<int> _endResidues;
        std::vector<std::vector<int>> _blockResidues;
        //std::vector<D_PResidues> _pResiduesBond;
        //std::vector<D_WResidues> _wResiduesBond;




        //std::vector<ModifiedAtomBondInfo> _atomsBondLoaded;
        int _bondCutoff = 3;// another case is 3: if bondCutoff is 3, the loop to find _exclusions runs twice to include particles that are 2 bonds away.
        std::vector<Coords3D> _totalForces;
        std::vector<Coords3D> _velocities;
        std::vector<double> _masses;
        std::vector<double> _inverseMasses;
        std::vector<double> _kineticEnergies;// a vector of _kineticEnergy for each atom
        //std::vector<double> _potentialEnergies;// a vector of _potentialEnergy for each atom does not exist since potential is calculated per force not per atom
        double _totalKEnergy;// total kinetic Energy of the system at current step
        double _totalPEnergy;// total potential Energy of the system at current step
        double _totalEnergy;// total Energy=potential+kinetic of the system at current step
        int _enabledForces; // Store enabled forces bitmask

        PeriodicBoundaryCondition::BoxInfo _boxInfo;


        std::vector<PDBAtom> _outputTemplate;//to store the exracted output.pdb format from input.pdb file 
        int _Modelnum;



        std::vector<double> _dt;
        double errorTol=0.0002;//check and adjust later 
        double maxStepSize=0.02;//check and adjust later this is in ps and is 20 fs


        // Force enable flags
        bool _harmonicBondForceEnabled;
        bool _harmonicAngleForceEnabled;
        bool _periodicTorsionForceEnabled;
        bool _nonbondedForceEnabled;

        //void InitializeFromFiles(const std::string& systemFilename, const std::string& stateFilename);

        //void Initialize(const vector<Coords3D>& atomPositions);
        //pair<vector<Coords3D>, vector<Coords3D>> Initialize(const vector<Coords3D>& atomPositions);
        //void Initialize(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
        //std::vector<Coords3D> CalculateForces(const std::vector<Coords3D>& atomPositions, const std::vector<PTorsionParams>& torsionParams, const std::vector<AngleParams>& angleParams, std::vector<Coords3D>& totalForces);
        

        // Allocate memory on CPU for conversion to double3
        double3* _atomPositions_double3 = nullptr;
        double3* _forces_double3 = nullptr;
        double3* _velocities_double3 = nullptr;
        double3* _boxSize_double3 = nullptr;
        double3* _lb_double3 = nullptr;
        double3* _ub_double3 = nullptr;
        double* _bondPEnergies = nullptr;



        // GPU data members
        double3* d_atomPositions;      // Positions of the atoms (GPU)
        double* d_masses;              // Masses of atoms (GPU)
        double* d_inverseMasses;              // Masses of atoms (GPU)
        double3* d_velocities;         // Velocities of atoms (GPU)
        double3* d_totalForces;        // Total forces on atoms (GPU)
        BondParams* d_bondParams;          // Bond parameters (GPU)
        PTorsionParams* d_torsionParams;       // Torsion parameters (GPU)
        AngleParams* d_angleParams;         // Angle parameters (GPU)
        NonbondedParams* d_nonbondedParams;     // Nonbonded interaction parameters (GPU)
        //ModifiedAtomBondInfo* d_atomsBondLoaded;
        D_Residues* d_residues = nullptr; //protein residue device pointer
        D_CudaBonds* d_cudaBonds = nullptr; //protein residue device pointer

        // CudaDatastructures
        //int* d_startResidues;
        //int* d_endResidues;

        double* d_totalPEnergy;
        double* d_bondPEnergies;
        double* d_kineticEnergies;
        double* d_totalKEnergy;
        double* d_totalEnergy;
        double3* d_boxsize;
        double3* d_lb;
        double3* d_ub;
        int* d_numAtoms;
        int* d_numBonds;
        int* d_numAngles;
        int* d_numTorsions;
        int* d_numNonbonded;
        
        
        
        
        
        
        
        void extractBoxBoundaries(const std::vector<Coords3D>& atomPositions, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        void CalculateForces();
        //pair<vector<Coords3D>, vector<Coords3D>> Integrate(int& StepNum, double& StepSize);
        //std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Integrate2(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& Step, double& StepSize);
        void Integrate(int& Step);
        void TotalEnergy(double& timestep);
        //void Report(const string& outputFilename, int step);
        void Report(const std::string& inputFilename, const std::string& outputFilename, int& step, double& timestep, int& interval);
        void CleanupGPU();

    };
}

#endif // CUDAENGINE_H
