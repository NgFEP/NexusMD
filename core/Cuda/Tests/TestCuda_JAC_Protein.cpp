#include <iostream>
#include <string>
#include "TaskDispatcher.h"

using namespace std;

int main() {
    // Simulation parameters
    string systemFilename = "system_JAC_Protein.xml";
    string stateFilename = "state_JAC_Protein.xml";
    string outputFilename = "output_SimRun_JAC_Protein_HBond.pdb";
    string inputFilename = "JAC_Protein.pdb";
    double stepSize = 0.001;//0.001 pico s =1 femto s
    int totalSteps = 10;
    int interval = 10;


    vector<Coords3D> atomPositions = {};

    vector<PTorsionParams> torsionParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force

    vector<BondParams> bondParams{};

    vector<AngleParams> angleParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force
    NonbondedParams nonbondedParams;// in this case nonbondedParams.particles = {} is automatically equal to zero

    vector<double> masses = {};  // Specific masses for the particles
    PeriodicBoundaryCondition::BoxInfo boxInfo;
    boxInfo.boxSize = {};


    // Create TaskDispatcher using the factory method
    auto dispatcher = TaskDispatcher::CreateDispatcher("CUDA");// Choose between "CPU" and "CUDA"

    try {

        //long startTime = clock();

        // Dispatch and run the simulation
        dispatcher->Simulate(systemFilename, stateFilename, inputFilename, outputFilename, stepSize, totalSteps, interval, atomPositions, masses, torsionParams, bondParams, angleParams, nonbondedParams, boxInfo);
        cout << "Simulation completed successfully." << endl;

        //long finishTime = clock();
        //cout << "Runtime took " << (finishTime - startTime) << " millis" << endl;

    }

    catch (const exception& e) {
        cerr << "Error during simulation: " << e.what() << endl;
    }

    return 0;
}






































////#include "stdafx.h"
//
//#pragma once
//
//#include <iostream>
//#include <vector>
//#include <string>
//#include <sstream>
//#include <fstream>
//#include <utility>
//#include "Engine.h"
//#include "Coords3D.h"
//
//using namespace std;
//using namespace BaseLine;
//
//int main() {
//    string systemFilename = "system_dialanine.xml"; // Path to the system file
//    string stateFilename = "state_dialanine.xml"; // Path to the state file
//    string outputFilename = "output_SimRun_Dialanine_PBC_test_PME_Bond_Angle.pdb"; // Output file for reporting
//    string inputFilename = "dialanine_VMD2.pdb"; // Output file for reporting
//
//    // Assuming these vectors need to be passed to the RunSimulation method.
//    //vector<Coords3D> totalForces; // Initialized empty, assuming Engine handles their setup.
//    //vector<Coords3D> velocities; // Initialized empty, similarly.
//
//    Engine engine(systemFilename, stateFilename, {}, {}, {}, {}, {});
//
//    double StepSize = 0.001; // Example step size
//    int TotalSteps = 5; // Example total number of steps for the simulation
//    int interval = 1; // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)
//    cout << "Total Frames #: " << abs(TotalSteps / interval) << endl;
//
//    // Run the simulation with specified parameters
//    try {
//        engine.RunSimulation(inputFilename, outputFilename, StepSize, TotalSteps, interval);// , systemFilename, stateFilename, totalForces, velocities);
//
//        cout << "Simulation completed successfully." << endl;
//    }
//    catch (const exception& e) {
//        cerr << "Error during simulation: " << e.what() << endl;
//    }
//
//    return 0;
//}
//
//










//#include "Engine.h"
//#include <iostream>
//
//using namespace BaseLine;
//
//int main() {
//
//    string systemFilename = "system.xml"; // Path to your system file
//    string stateFilename = "state.xml"; // Path to your state file
//    string outputFilename = "output.txt"; // Output file for reporting
//    Engine engine(systemFilename, stateFilename);
//
//    //Engine engine;
//
//    double StepSize = 0.001; // Example step size
//    int TotalSteps = 5; // Example total number of steps for the simulation
//    int Step = 0; // Starting step number, assuming this is intended to be an iteration count or similar
//
//    // Run the simulation with specified parameters
//    try {
//        engine.RunSimulation(outputFilename, StepSize, TotalSteps, systemFilename, stateFilename, );
//
//        cout << "Simulation completed successfully." << endl;
//    }
//    catch (const exception& e) {
//        cerr << "Error during simulation: " << e.what() << endl;
//    }
//
//    return 0;
//}