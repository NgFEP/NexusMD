#ifndef CUDAFORCES_H
#define CUDAFORCES_H
//This class aimes to calculate totalforces on all the atoms by going through and calling various forces parameters inside the system.xml

#include "Coords3D.h"
#include <set>
#include "CudaStateXMLParser.h"
#include "SystemXMLParser.h"
#include "CudaPeriodicTorsionForce.h"
#include "CudaHarmonicAngleForce.h"
#include "CudaNonbondedForce.h"
#include "PeriodicBoundaryCondition.h"
#include "HarmonicBondForceKernel.h"

//#include <string>
#include <vector>

namespace Cuda {

    class Forces {
    public:
        //adding PeriodicTorsionForce to the total forces by going through each torsion and calculates each force on involved atoms and adding them to the total forces vector
        static void AddPTorsion(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<PTorsionParams>& torsionParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        static void AddHBond(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<BondParams>& bondParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo); // New function declaration
        static void AddHAngle(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<AngleParams>& angleParams, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo);
        static void AddNonBondElectroPME(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const NonbondedParams& params, double& totalPEnergy, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const std::vector<std::set<int>>& exclusions);

    private:


    };

} // namespace Cuda

#endif // CUDAFORCES_H
