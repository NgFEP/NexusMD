#include <iostream>
#include <string>
#include "TaskDispatcher.h"

using namespace std;

int main() {
    // Simulation parameters
    string systemFilename = "system_waterbox.xml";
    string stateFilename = "state_waterbox.xml";
    string outputFilename = "output_SimRun_waterbox.pdb";
    string inputFilename = "waterbox.pdb";
    double stepSize = 0.001; // 0.001 pico s = 1 femto s
    int totalSteps = 2000;
    int interval = 100;

    // Create TaskDispatcher using the factory method for CPU
    auto dispatcher = TaskDispatcher::CreateDispatcher("cpu");

    // Directly enable forces
    dispatcher->EnableHarmonicBondForce();
    //dispatcher->EnableHarmonicAngleForce();
    //dispatcher->EnablePeriodicTorsionForce();
    //dispatcher->EnableNonbondedForce();
    //Thermostat
    //dispatcher->AndersenThermostat(275,11);//CollisionFrequency, Temperature

    dispatcher->AndersenThermostat(300,10);//CollisionFrequency, Temperature
    //dispatcher->MonteCarloBarostat(1.0, 25);

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
