#ifndef NONBONDEDFORCE_H
#define NONBONDEDFORCE_H

#include "Coords3D.h"
#include "Types.h"
#include <vector>
#include <set>
#include "SystemXMLParser.h" // This includes the definition for NonbondedParams
#include "PeriodicBoundaryCondition.h"
#include "NeighborList.h"


#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef SQRT_PI
#define SQRT_PI 1.77245385090551602729
#endif

#ifndef TWO_OVER_SQRT_PI
#define TWO_OVER_SQRT_PI 1.1283791670955126 // Accurate value of 2/sqrt(pi)
#endif

//#define ONE_4PI_EPS0 8.9875517873681764e9 //​ unit N.m2/c2
#ifndef ONE_4PI_EPS0
#define ONE_4PI_EPS0 138.93545764446428693393 // kJ⋅nm/(mol⋅e2)
#endif







namespace BaseLine {


    class NonbondedForce {
    public:
        // Constructor
        NonbondedForce();


        // To calculate forces for all particles based on their positions and nonbonded parameters
        std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const std::vector<std::set<int>>& exclusions);



        // To apply exclusions to adjust forces and potential energy calculations
        //static void applyExclusions(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& forces, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxBuffer);

    private:
        Coords3D recipBoxVec;
        void initializeReciprocalVectors(const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        std::vector<int> gridSize = { 0, 0, 0 };
        int findFFTDimension(int minimum);
        void initializeGridSize(const PeriodicBoundaryCondition::BoxInfo& boxInfo);


        std::vector<std::vector<double>> bsplines_moduli;
        std::vector<std::vector<std::vector<double>>> bsplines_theta;
        std::vector<std::vector<std::vector<double>>> bsplines_dtheta;

        std::vector<std::vector<int>> particleIndex;
        std::vector<std::vector<double>> particleFraction;
        double epsilon_r = 1.0;    //epsilon_r   Dielectric coefficient, typically 1.0.
        NonbondedParams _ModifiedNonbondedParams; // used for direct space calculations

        
                void InitparticleIndexAndFraction(const std::vector<Coords3D>& atomPositions);
        void computeBSplineParameters(int pmeOrder, const std::vector<Coords3D>& atomPositions);


        double alphaEwald;  // Ewald parameter, adjust as necessary
        const int Num_Table_Points = 2048;  // Number of points in the table
        std::vector<double> erfcTable;
        double RECIP_EXP_FACTOR;


        void initializeErfcTable(const NonbondedParams& params);

        // To apply periodic boundary conditions to coordinates
        void applyPeriodicBoundary(Coords3D& pos, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        void findAtomGridIndex(const std::vector<Coords3D>& atomPositions, std::vector<Int3D>& gridIndices, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo);


        // To spread charges to the PME grid
        void gridSpreadCharge_Reference(const std::vector<Coords3D>& atomPositions, ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const int& pmeOrder);
        void perform3DFFT(ComplexGrid& pmeGrid, bool forward);



        void InitializeAlphaEwald(const NonbondedParams& params);

        void reciprocalConvolutionAndEnergy(ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& energy);

        void PMEGridForces(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy, const int& pmeOrder);
        void ReciprocalPMEcalculateForcesAndEnergy(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        void NonbondedParamsModifier(const std::vector<Coords3D>& atomPositions, const NonbondedParams& params);
        void CalculateSelfEnergy(const NonbondedParams& params, double& totalPEnergy);
        void DirectForcesAndEnergy(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const std::vector<std::set<int>>& exclusions, double& totalPEnergy);
        void CalculateExclusionEnergyAndForces(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const std::vector < std::set<int >> &exclusions, double& totalPEnergy);

    };

} // namespace BaseLine

#endif // NONBONDEDFORCE_H
