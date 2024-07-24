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
#include "KineticEnergy.h"
#include "PeriodicBoundaryCondition.h"
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include "Coords3D.h"
#include "Reporter.h"

namespace BaseLine {

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
            const std::vector<Coords3D>& atomPositions = std::vector<Coords3D>(),
            const std::vector<double>& masses = std::vector<double>(),
            const std::vector<PTorsionParams>& torsionParams = std::vector<PTorsionParams>(),
            const std::vector<HBondParams>& bondParams = std::vector<HBondParams>(),
            const std::vector<HAngleParams>& angleParams = std::vector<HAngleParams>(),
            const NonbondedParams& nonbondedParams = NonbondedParams(),
            const PeriodicBoundaryCondition::BoxInfo& boxInfo = PeriodicBoundaryCondition::BoxInfo() // Default empty box info
        );
        
        
        //void Inputs(const string& systemFilename, const string& stateFilename);
        //void RunSimulation(const string& outputFilename, double timestep, int numSteps, const string& systemFilename, const string& stateFilename);
        void RunSimulation(const std::string& inputFilename, const std::string& outputFilename, double& timestep, int& numSteps, int& interval);// , const std::string& systemFilename, const std::string& stateFilename, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);


    private:
        void InitializeSimulationParameters();
        bool _RealSimRun=false;
        std::size_t _numAtoms;
        std::string _systemFilename;
        std::string _stateFilename;
        std::vector<Coords3D> _atomPositions;
        std::vector<PTorsionParams> _torsionParams;
        std::vector<HBondParams> _bondParams;
        std::vector<HAngleParams> _angleParams;
        NonbondedParams _nonbondedParams;
        std::vector<Coords3D> _totalForces;
        std::vector<Coords3D> _velocities;
        std::vector<double> _masses;
        std::vector<double> _inverseMasses;
        std::vector<double> _kineticEnergies;// a vector of _kineticEnergy for each atom
        //std::vector<double> _potentialEnergies;// a vector of _potentialEnergy for each atom does not exist since potential is calculated per force not per atom
        double _totalKEnergy;// total kinetic Energy of the system at current step
        double _totalPEnergy;// total potential Energy of the system at current step
        double _totalEnergy;// total Energy=potential+kinetic of the system at current step

        PeriodicBoundaryCondition::BoxInfo _boxInfo;


        std::vector<PDBAtom> _outputTemplate;//to store the exracted output.pdb format from input.pdb file 
        std::size_t _Modelnum;



        std::vector<double> _dt;
        double errorTol=0.0002;//check and adjust later 
        double maxStepSize=0.02;//check and adjust later this is in ps and is 20 fs

        //void InitializeFromFiles(const std::string& systemFilename, const std::string& stateFilename);

        //void Initialize(const vector<Coords3D>& atomPositions);
        //pair<vector<Coords3D>, vector<Coords3D>> Initialize(const vector<Coords3D>& atomPositions);
        //void Initialize(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);
        //std::vector<Coords3D> CalculateForces(const std::vector<Coords3D>& atomPositions, const std::vector<PTorsionParams>& torsionParams, const std::vector<HAngleParams>& angleParams, std::vector<Coords3D>& totalForces);
        void extractBoxBoundaries(const std::vector<Coords3D>& atomPositions, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        void CalculateForces();
        //pair<vector<Coords3D>, vector<Coords3D>> Integrate(int& StepNum, double& StepSize);
        //std::pair<std::vector<Coords3D>, std::vector<Coords3D>> Integrate2(std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& velocities, std::vector<Coords3D>& totalForces, std::vector<double>& masses, int& Step, double& StepSize);
        void Integrate(int& Step);
        void TotalEnergy(double& timestep);
        //void Report(const string& outputFilename, int step);
        void Report(const std::string& inputFilename, const std::string& outputFilename, int& step, int& interval);

    };
}

#endif // ENGINE_H
