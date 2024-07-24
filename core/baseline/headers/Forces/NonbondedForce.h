#ifndef NONBONDEDFORCE_H
#define NONBONDEDFORCE_H

#include "Coords3D.h"
#include "Types.h"
#include <vector>
#include "SystemXMLParser.h" // This includes the definition for NonbondedParams
#include "PeriodicBoundaryCondition.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace BaseLine {


    //struct NeighborList {
    //    std::vector<pair<int, int>> pairs;
    //};

    class NonbondedForce {
    public:
        // Constructor
        NonbondedForce();

        // Calculate forces for all particles based on their positions and nonbonded parameters
        std::vector<Coords3D> calculateForces(const std::vector<Coords3D>& atomPositions, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        //void DirectAndReciprocalForcesAndEnergy(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy, const int& pmeOrder);
        void ReciprocalForcesAndEnergy(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy, const int& pmeOrder);
        //void DirectForcesAndEnergy(std::vector<Coords3D>& forces, const std::vector<Coords3D>& atomPositions, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy);



        // Apply exclusions to adjust forces and potential energy calculations
        //static void applyExclusions(const std::vector<Coords3D>& atomPositions, std::vector<Coords3D>& forces, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxBuffer);

    private:
        //double recipBoxVecX, recipBoxVecY, recipBoxVecZ;
        Coords3D recipBoxVec;
        //int pmeOrder;
        void initializeReciprocalVectors(const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        std::vector<int> gridSize = { 0, 0, 0 };
        int findFFTDimension(int minimum);
        void initializeGridSize(const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        //std::vector<double> bSplineModuliX;
        //std::vector<double> bSplineModuliY;
        //std::vector<double> bSplineModuliZ;

        std::vector<std::vector<double>> bsplines_moduli;
        std::vector<std::vector<std::vector<double>>> bsplines_theta;
        std::vector<std::vector<std::vector<double>>> bsplines_dtheta;

        std::vector<std::vector<int>> particleIndex;
        std::vector<std::vector<double>> particleFraction;
        double epsilon_r = 1.0;    //epsilon_r   Dielectric coefficient, typically 1.0.

        

        //std::vector<double> computeBSplineModuli(int pmeOrder);
        void InitparticleIndexAndFraction(const std::vector<Coords3D>& atomPositions);
        void computeBSplineParameters(int pmeOrder, const std::vector<Coords3D>& atomPositions);

        //void initializeBSplineModuli(int pmeOrder);

        double alphaEwald;  // Ewald parameter, adjust as necessary
        const int Num_Table_Points = 2048;  // Number of points in the table
        std::vector<double> erfcTable;
        double RECIP_EXP_FACTOR;


        void initializeErfcTable(const NonbondedParams& params);

        // Apply periodic boundary conditions to coordinates
        void applyPeriodicBoundary(Coords3D& pos, const PeriodicBoundaryCondition::BoxInfo& boxInfo);

        void findAtomGridIndex(const std::vector<Coords3D>& atomPositions, std::vector<Int3D>& gridIndices, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo);


        // Spread charges to the PME grid
        void gridSpreadCharge_Reference(const std::vector<Coords3D>& atomPositions, ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const int& pmeOrder);
        void perform3DFFT(ComplexGrid& pmeGrid, bool forward);


        //void finishSpreadCharge(ComplexGrid& pmeGrid, const int& pmeOrder);

        void InitializeAlphaEwald(const NonbondedParams& params);

        void reciprocalConvolution_reference(ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& energy);

        // void computeNeighborList(NeighborList& neighborList, const std::vector<Coords3D>& atomPositions, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double maxDistance, double minDistance, bool usePeriodic);

    
    };

} // namespace BaseLine

#endif // NONBONDEDFORCE_H
