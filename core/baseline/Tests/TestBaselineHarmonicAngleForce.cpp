//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "Engine.h"
#include "Coords3D.h"

using namespace std;
using namespace BaseLine;

int main() {

    // Assuming these vectors need to be passed to the RunSimulation method.
    //vector<Coords3D> totalForces; // Initialized empty, assuming Engine handles their setup.
    //vector<Coords3D> velocities; // Initialized empty, similarly.


    // first harmonic angle force info: 
    // <Angle a="1.911135530933791" k="418.40000000000003" p1="0" p2="4" p3="5"/>

    vector<Coords3D> atomPositions = {
        {2.5160000324249268, 1.4160000085830688, 1.9440000057220459},
        {2.509000062942505, 1.3919999599456787, 1.7979999780654907},
        {2.4700000286102295, 1.2929999828338623, 1.7760000228881836}

    };

    vector<PTorsionParams> torsionParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic angle force
    vector<HBondParams> bondParams{};//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant
    vector<HAngleParams> angleParams{ {0, 1, 2, 1.911135530933791, 418.40000000000003} };//p1,p2,p3,a,k //double a; // Ideal angle in radians and double k; // Force constant
    vector<double> masses = { 14.01, 12.01, 1.008 };
    
    Engine engine("", "", atomPositions, masses, torsionParams, bondParams, angleParams);//engine construction for test run
    string outputFilename = "output_test_HAF_1Angle.txt"; // Output file for reporting

    double StepSize = 0.001; // Example step size
    int TotalSteps = 3; // Example total number of steps for the simulation
    int interval = 1;

    // Run the simulation with specified parameters
    try {
        //engine.RunSimulation(outputFilename, StepSize, TotalSteps);// , systemFilename, stateFilename, totalForces, velocities);
        engine.RunSimulation("", outputFilename, StepSize, TotalSteps, interval);// , systemFilename, stateFilename, totalForces, velocities);

        cout << "Simulation completed successfully." << endl;
    }
    catch (const exception& e) {
        cerr << "Error during simulation: " << e.what() << endl;
    }

    return 0;
}