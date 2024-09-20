//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "Coords3D.h"
#include <cstring>
#include <cstdio>
#include "StateXMLParser.h"
#include "VerletIntegration.h"
#include "PeriodicBoundaryCondition.h"
#include <cmath> // Include for mathematical functions like sqrt()
#include <algorithm>


using namespace BaseLine;
using namespace std;



// Define the static member variable with namespace qualification
vector<Coords3D> VerletIntegration::UpdatedAtomPositions;
VerletIntegration::VerletIntegration() {
    // Constructor implementation
}

VerletIntegration::~VerletIntegration() {
    // Destructor implementation
}


void VerletIntegration::InverseMasses(vector<double>& masses, vector<double>& inverseMasses) {
    // Initialize inverse masses only once to prevent overcomputation of inverseMasses for every steps
    if (inverseMasses.empty()) {
        inverseMasses.resize(masses.size());
        for (int i = 0; i < masses.size(); ++i) {
            inverseMasses[i] = (masses[i] == 0.0) ? 0.0 : 1.0 / masses[i];
        }
        //cout << "inverseMasses: " << inverseMasses[0] << "   " << inverseMasses[1] << "   " << inverseMasses[2] << "   ";
    }
}


//modified version NexaBind modified with cuda openmm method
void VerletIntegration::Advance(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& inverseMasses, int& StepNum, const vector<double>& dt, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    double dtPos = dt[1];
    double dtVel = 0.5 * (dt[0] + dt[1]);
    double scale = dtVel;// / (double)0x100000000;
    
    int numberOfAtoms = atomPositions.size();

    UpdatedAtomPositions.assign(numberOfAtoms, Coords3D(0, 0, 0));

  //vector<Coords3D>& velocities2= velocities;

    //if (StepNum==1) {
    //    cout << "";
    //}
    // Perform the integration
    for (int i = 0; i < numberOfAtoms; ++i) {
        if (inverseMasses[i] != 0.0) {
            for (int j = 0; j < 3; ++j) {
                //note that the following calculated velocity is at t+StepSize/2 and not t+StepSize
                velocities[i][j] += inverseMasses[i] * totalForces[i][j] * scale;

                if (velocities[i][j] > 200) {
                    cout << "";
                }

                //note that the following calculated atompositions is at t+StepSize
                UpdatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j] * dtPos;

                // For PBC: Wrap coordinates to stay within box dimensions
                if (UpdatedAtomPositions[i][j] < boxInfo.lb[j]) {
                    UpdatedAtomPositions[i][j] += boxInfo.boxSize[j];
                }
                else if (UpdatedAtomPositions[i][j] > boxInfo.ub[j]) {
                    UpdatedAtomPositions[i][j] -= boxInfo.boxSize[j];
                }
                // loop version which is not desired
                //while (true) {
                    //if (UpdatedAtomPositions[i][j] < boxInfo.lb[j]) {
                    //    UpdatedAtomPositions[i][j] += boxInfo.boxSize[j];
                    //}
                    //else if (UpdatedAtomPositions[i][j] > boxInfo.ub[j]) {
                    //    UpdatedAtomPositions[i][j] -= boxInfo.boxSize[j];
                    //}
                    // loop until the updated position is located in the sim box 
                    //if (UpdatedAtomPositions[i][j] >= boxInfo.lb[j] && UpdatedAtomPositions[i][j] <= boxInfo.ub[j]) {
                    //    break;
                    //}
                //}

                // ***create an error to ask the user for smaller step sizes if particles move more than one box size at one step***



            }
        }
    }
    double oneOverDt = 1.0 / dt[1];

    for (int i = 0; i < numberOfAtoms; ++i) {
        if (inverseMasses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {
                // now these velocities and atomPositions are at the next step of t+StepSize 
                

                // bug fixed: due to PBC updatedAtomPositions is updated and atompositions is still on the other side of the box which creats very large velocity values, since the following equation is the same as                 UpdatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j] * dtPos;
                // therfore there is no point in adding it.
                //velocities[i][j] = (UpdatedAtomPositions[i][j] - atomPositions[i][j]) * oneOverDt;
                
                //if (velocities[i][j] > 200) {
                //    cout << "";
                //}
                atomPositions[i][j] = UpdatedAtomPositions[i][j];
            }
    }

}



////openmm cuda version
//void VerletIntegration::Advance(const vector<double>& dt, vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& inverseMasses, vector<Coords3D>& posDelta) {
//    double dtPos = dt[1];
//    double dtVel = 0.5 * (dt[0] + dt[1]);
//    double scale = dtVel / (double)0x100000000;
//
//
//    for (size_t i = 0; i < atomPositions.size(); ++i) {
//        if (inverseMasses[i] != 0.0) {
//            velocities[i][0] += scale * totalForces[i][0] * inverseMasses[i];
//            velocities[i][1] += scale * totalForces[i][1] * inverseMasses[i];
//            velocities[i][2] += scale * totalForces[i][2] * inverseMasses[i];
//            posDelta[i][0] = velocities[i][0] * dtPos;
//            posDelta[i][1] = velocities[i][1] * dtPos;
//            posDelta[i][2] = velocities[i][2] * dtPos;
//            atomPositions[i][0] += posDelta[i][0];
//            atomPositions[i][1] += posDelta[i][1];
//            atomPositions[i][2] += posDelta[i][2];
//        }
//    }
//
//
//    double oneOverDt = 1.0 / dt[1];
//
//    for (size_t i = 0; i < atomPositions.size(); ++i) {
//        if (inverseMasses[i] != 0.0) {
//            velocities[i][0] = posDelta[i][0] * oneOverDt;
//            velocities[i][1] = posDelta[i][1] * oneOverDt;
//            velocities[i][2] = posDelta[i][2] * oneOverDt;
//        }
//    }
//
//    //// Perform the integration
//    //for (int i = 0; i < numberOfAtoms; ++i) {
//    //    if (masses[i] != 0.0) {
//    //        for (int j = 0; j < 3; ++j) {
//    //            //note that the following calculated velocity is at t+StepSize/2 and not t+StepSize
//    //            velocities[i][j] += inverseMasses[i] * totalForces[i][j] * StepSize / 2.0;
//    //            //note that the following calculated atompositions is at t+StepSize
//    //            UpdatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j] * StepSize;
//    //        }
//    //    }
//    //}
//
//
//}

//void VerletIntegration::integrateVerletPart2(const vector<double>& dt, vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<double>& inverseMasses, const vector<Coords3D>& posDelta) {
//    double oneOverDt = 1.0 / dt[1];
//
//    for (size_t i = 0; i < atomPositions.size(); ++i) {
//        if (inverseMasses[i] != 0.0) {
//            velocities[i][0] = posDelta[i][0] * oneOverDt;
//            velocities[i][1] = posDelta[i][1] * oneOverDt;
//            velocities[i][2] = posDelta[i][2] * oneOverDt;
//        }
//    }
//
//    //for (int i = 0; i < numberOfAtoms; ++i) {
//    //    if (masses[i] != 0.0)
//    //        for (int j = 0; j < 3; ++j) {
//    //            //now these velocities and atomPositions are at the next step of t+StepSize 
//    //            velocities[i][j] = (UpdatedAtomPositions[i][j] - atomPositions[i][j]) / StepSize;
//    //            atomPositions[i][j] = UpdatedAtomPositions[i][j];
//    //        }
//    //}
//
//
//
//
//}
// dynamically adjust the integration step size based on the forces experienced by the particles and their masses
void VerletIntegration::SelectVerletStepSize(vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& inverseMasses, vector<double>& dt, double errorTol, double maxStepSize) {
    double totalError = 0.0;
    double scale = 1.0; // / (double)0x100000000;

    for (int i = 0; i < velocities.size(); ++i) {
        double invMass = inverseMasses[i];
        if (invMass != 0.0) {
            double err = (scale * totalForces[i][0]) * (scale * totalForces[i][0]) +
                (scale * totalForces[i][1]) * (scale * totalForces[i][1]) +
                (scale * totalForces[i][2]) * (scale * totalForces[i][2]);
            err *= invMass * invMass;
            totalError += err;
        }
    }
    
    totalError = sqrt(totalError / (3 * velocities.size()));
    double newStepSize = sqrt(errorTol / totalError);
    double oldStepSize = dt[1];
    if (oldStepSize > 0) {
        newStepSize = min(newStepSize, oldStepSize * 2.0);
    }
    if (newStepSize > oldStepSize && newStepSize < 1.1 * oldStepSize) {
        newStepSize = oldStepSize;
    }
    if (newStepSize > maxStepSize) {
        newStepSize = maxStepSize;
    }
    dt[1] = newStepSize;
}





//// Define the static member variable with namespace qualification
//vector<double> VerletIntegration::inverseMasses;
//vector<Coords3D> VerletIntegration::updatedAtomPositions;
//
//VerletIntegration::VerletIntegration() {}
//
//VerletIntegration::~VerletIntegration() {}
//
//void VerletIntegration::advance(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces, vector<double>& masses, double stepSize) {
//    int numberOfAtoms = atomPositions.size();
//    updatedAtomPositions.assign(numberOfAtoms, Coords3D(0, 0, 0));
//
//    // Initialize inverse masses only once at the beginning of the simulation
//    if (inverseMasses.empty()) {
//        inverseMasses.resize(masses.size());
//        for (size_t i = 0; i < masses.size(); ++i) {
//            inverseMasses[i] = (masses[i] == 0.0) ? 0.0 : 1.0 / masses[i];
//        }
//    }
//
//    integrateVerletPart1(atomPositions, velocities, totalForces);
//    integrateVerletPart2(atomPositions, velocities);
//}
//
//void VerletIntegration::integrateVerletPart1(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities, vector<Coords3D>& totalForces) {
//    double dtVel = 0.5;  // Placeholder for delta time for velocity calculation
//    for (int i = 0; i < atomPositions.size(); ++i) {
//        if (inverseMasses[i] != 0.0) {
//            for (int j = 0; j < 3; ++j) {
//                velocities[i][j] += inverseMasses[i] * totalForces[i][j] * dtVel;
//                updatedAtomPositions[i][j] = atomPositions[i][j] + velocities[i][j];
//            }
//        }
//    }
//}
//
//void VerletIntegration::integrateVerletPart2(vector<Coords3D>& atomPositions, vector<Coords3D>& velocities) {
//    double dtPos = 1.0;  // Placeholder for delta time for position update
//    for (int i = 0; i < atomPositions.size(); ++i) {
//        if (inverseMasses[i] != 0.0) {
//            for (int j = 0; j < 3; ++j) {
//                velocities[i][j] = (updatedAtomPositions[i][j] - atomPositions[i][j]) / dtPos;
//                atomPositions[i][j] = updatedAtomPositions[i][j];
//            }
//        }
//    }
//}
//
