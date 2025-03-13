#ifndef ENSEMBLEPARAMETERS_H
#define ENSEMBLEPARAMETERS_H

#include <vector>
#include "Coords3D.h"  // Assuming Coords3D is defined elsewhere in your code
#include <iostream>
#include "DataStructures.h"


#define boltz (0.008314462618)	/*kJ / (mol·K)*/
//#define AMU_TO_GRAM (1.66054e-24);  // Conversion factor: 1 amu = 1.66054 × 10-24 g
#define NM3_TO_ML (1e-21);
#define AVOGADRO (6.02214076e23)          

namespace BaseLine {


    // Class to compute temperature based on velocities and kinetic energy
    class EnsembleParameters{
    public:
        EnsembleParameters();
        ~EnsembleParameters();
    
        // Updated ComputeTemperature to take kineticEnergy as input
        double ComputeTemp(int& numAtoms, double& kineticEnergy, bool& removeCOM, const std::vector<Constraint>& constraints);
        double ComputeTotalMass(const std::vector<double>& Masses);
        double ComputeVolume(const Coords3D& boxSize);
        double ComputeDensity(double& TotalMass, double& Volume);

    private:
        // Function to calculate DOF for the system based on constraints
        int numDOF(int& numAtoms, bool& removeCOM, const std::vector<Constraint>& constraints);
    };

}

#endif // ENSEMBLEPARAMETERS_H
