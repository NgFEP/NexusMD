//Engine is the heart of the Nexus which does the actual simulation (evolving the system by combining the OpenMM system data, accumulating forces, and verlet integrator)
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
#include <set>
#include "Coords3D.h"
#include "Reporter.h"
#include "Exclusions.h"


namespace BaseLine {




    class Engine {
    public:
        // Constructor declaration
        Engine(//this decleration allows optional inputs for engine object construction useful for both test and production runs
            const std::string& systemFilename = "",
            const std::string& stateFilename = "",
            const std::vector<Coords3D>& atomPositions = std::vector<Coords3D>(),
            const std::vector<double>& masses = std::vector<double>(),
            const std::vector<PTorsionParams>& torsionParams = std::vector<PTorsionParams>(),
            const std::vector<BondParams>& bondParams = std::vector<BondParams>(),
            const std::vector<AngleParams>& angleParams = std::vector<AngleParams>(),
            const NonbondedParams& nonbondedParams = NonbondedParams(),
            const PeriodicBoundaryCondition::BoxInfo& boxInfo = PeriodicBoundaryCondition::BoxInfo() // Default empty box info
        );
        
        
        void RunSimulation(const std::string& inputFilename, const std::string& outputFilename, double& timestep, int& numSteps, int& interval);// , const std::string& systemFilename, const std::string& stateFilename, std::vector<Coords3D>& totalForces, std::vector<Coords3D>& velocities);


    private:
        void InitializeSimulationParameters();
        bool _RealSimRun=false;
        int _numAtoms;
        std::string _systemFilename;
        std::string _stateFilename;
        std::vector<Coords3D> _atomPositions;
        std::vector<PTorsionParams> _torsionParams;
        std::vector<BondParams> _bondParams;
        std::vector<AngleParams> _angleParams;
        NonbondedParams _nonbondedParams;
        std::vector<std::set<int>> _exclusions;
        int _bondCutoff = 3;// if bondCutoff is 3, the loop to find _exclusions runs twice to include particles that are 2 bonds away.
        std::vector<Coords3D> _totalForces;
        std::vector<Coords3D> _velocities;
        std::vector<double> _masses;
        std::vector<double> _inverseMasses;
        std::vector<double> _kineticEnergies;// a vector of _kineticEnergy for each atom
        double _totalKEnergy;// total kinetic Energy of the system at current step
        double _totalPEnergy;// total potential Energy of the system at current step
        double _totalEnergy;// total Energy=potential+kinetic of the system at current step

        PeriodicBoundaryCondition::BoxInfo _boxInfo;


        std::vector<PDBAtom> _outputTemplate;//to store the exracted output.pdb format from input.pdb file 
        std::size_t _Modelnum;



        std::vector<double> _dt;
        double errorTol=0.0002;//check and adjust later 
        double maxStepSize=0.02;//in ps = 20 fs

        void CalculateForces();
        void Integrate(int& Step);
        void TotalEnergy(double& timestep);
        void Report(const std::string& inputFilename, const std::string& outputFilename, int& step, double& timestep, int& interval);

    };
}

#endif // ENGINE_H
