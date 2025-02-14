#include <iostream>
#include <string>
#include "TaskDispatcher.h"

using namespace std;

int main() {
    // Simulation parameters
    string systemFilename = "system_Sample_Protein.xml";
    string stateFilename = "state_Sample_Protein.xml";
    string outputFilename = "output_SimRun_Sample_Protein.pdb";
    string inputFilename = "Sample_Protein.pdb";
    double stepSize = 0.001; // 0.001 pico s = 1 femto s
    int totalSteps = 10;
    int interval = 10;

    // Create TaskDispatcher using the factory method for CPU
    auto dispatcher = TaskDispatcher::CreateDispatcher("CUDA");

    // Directly enable forces
    dispatcher->EnableHarmonicBondForce();
    //dispatcher->EnableHarmonicAngleForce();
    //dispatcher->EnablePeriodicTorsionForce();
    //dispatcher->EnableNonbondedForce();

    try {
        // Dispatch and run the simulation
        dispatcher->Simulate(systemFilename, stateFilename, inputFilename, outputFilename,
            stepSize, totalSteps, interval);

        cout << "Simulation completed successfully with all forces enabled." << endl;
    }
    catch (const exception& e) {
        cerr << "Error during simulation: " << e.what() << endl;
    }

    return 0;
}
