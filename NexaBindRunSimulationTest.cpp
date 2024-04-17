//#include "stdafx.h"
//#include "Engine.h"
//#include "Coords3D.h"
//
//using namespace std;
//using namespace BaseLine;
//
//int main() {
//    string systemFilename = "system.xml"; // Path to the system file
//    string stateFilename = "state.xml"; // Path to the state file
//    string outputFilename = "output.txt"; // Output file for reporting
//
//    // Assuming these vectors need to be passed to the RunSimulation method.
//    vector<Coords3D> totalForces; // Initialized empty, assuming Engine handles their setup.
//    vector<Coords3D> velocities; // Initialized empty, similarly.
//
//    Engine engine(systemFilename, stateFilename);
//
//    double StepSize = 0.001; // Example step size
//    int TotalSteps = 2; // Example total number of steps for the simulation
//
//    // Run the simulation with specified parameters
//    try {
//        engine.RunSimulation(outputFilename, StepSize, TotalSteps);// , systemFilename, stateFilename, totalForces, velocities);
//
//        cout << "Simulation completed successfully." << endl;
//    }
//    catch (const exception& e) {
//        cerr << "Error during simulation: " << e.what() << endl;
//    }
//
//    return 0;
//}































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