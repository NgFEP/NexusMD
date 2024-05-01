#ifndef FORCES_H
#define FORCES_H
//This class aimes to calculate totalforces on all the atoms by going through and calling various forces parameters inside the system.xml

#include "Coords3D.h"
#include "StateXMLParser.h"
#include "SystemXMLParser.h"
#include "PeriodicTorsionForce.h"
#include "HarmonicBondForce.h" // To use the HarmonicBondForce class
#include "HarmonicAngleForce.h"
//#include <string>
#include <vector>

namespace BaseLine {

    class Forces {
    public:
        //adding PeriodicTorsionForce to the total forces by going through each torsion and calculates each force on involved atoms and adding them to the total forces vector
        static void AddPTorsion(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<PTorsionParams>& torsionParams, double& totalPEnergy);
        static void AddHBond(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<HBondParams>& bondParams, double& totalPEnergy); // New function declaration
        static void AddHAngle(std::vector<Coords3D>& totalForces, const std::vector<Coords3D>& atomPositions, const std::vector<HAngleParams>& angleParams, double& totalPEnergy);

    private:


    };

} // namespace BaseLine

#endif // FORCES_H
