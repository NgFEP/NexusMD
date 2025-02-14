//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "TaskDispatcher.h"
#include "SystemXMLParser.h"
#include "PeriodicBoundaryCondition.h"

using namespace std;

int main() {

    // Assuming these vectors need to be passed to the RunSimulation method.
    //vector<Coords3D> totalForces; // Initialized empty, assuming Engine handles their setup.
    //vector<Coords3D> velocities; // Initialized empty, similarly.


    // first harmonic bond force info: 
    // <Bond d=".1522" k="265265.6" p1="19" p2="4"/>

    vector<Coords3D> atomPositions = {
        {2.6519999504089355, 1.3969999551773071, 1.7430000305175781},//		<Position p1=19/>
        {2.509000062942505, 1.3919999599456787, 1.7979999780654907}//		<Position p2=4/>
    };

    vector<PTorsionParams> torsionParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force
    //vector<HBondParams> bondParams{ {0, 1, .1522000000000000, 265265.6000000000} };//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant
    vector<HBondParams> bondParams{ {0, 1, 0.1521999984979630, 265265.5937500000000000} };//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant 
    vector<HAngleParams> angleParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force
    NonbondedParams nonbondedParams;// in this case nonbondedParams.particles = {} is automatically equal to zero

    vector<double> masses = { 12.01000000000000, 12.01000000000000 };

    string inputFilename = "";
    string outputFilename = "output_test_HBF_1Bond.txt"; // Output file for reporting
    string systemFilename = "";
    string stateFilename = "";

    double stepSize = 0.001; // Example step size
    int totalSteps = 1; // Example total number of steps for the simulation
    int interval = 1;

    PeriodicBoundaryCondition::BoxInfo boxInfo;
    boxInfo.boxSize = { 3.0, 3.0, 3.0 }; // Example manual setting // note that 
    //// Run the simulation with specified parameters
    //try {
    //    //engine.RunSimulation(outputFilename, StepSize, TotalSteps);// , systemFilename, stateFilename, totalForces, velocities);
    //    engine.RunSimulation("", outputFilename, StepSize, TotalSteps, interval);// , systemFilename, stateFilename, totalForces, velocities);

    //    cout << "Simulation completed successfully." << endl;
    //}


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





//
//#include <iostream>
//#include <string>
//#include "TaskDispatcher.h"
//
//using namespace std;
//
//int main() {
//    // Simulation parameters
//    string systemFilename = "system_dialanine.xml";
//    string stateFilename = "state_dialanine.xml";
//    string outputFilename = "output_SimRun_Dialanine_PBC_test_PME_Bond_Angle.pdb";
//    string inputFilename = "dialanine_VMD2.pdb";
//    double stepSize = 0.001;
//    int totalSteps = 5;
//    int interval = 1;
//
//    // Create TaskDispatcher using the factory method
//    auto dispatcher = TaskDispatcher::CreateDispatcher("CUDA");// Choose between "CPU" and "CUDA"
//
//    try {
//
//        long startTime = clock();
//
//        // Dispatch and run the simulation
//        dispatcher->Simulate(systemFilename, stateFilename, inputFilename, outputFilename, stepSize, totalSteps, interval);
//        cout << "Simulation completed successfully." << endl;
//
//        long finishTime = clock();
//        cout << "Runtime took " << (finishTime - startTime) << " millis" << endl;
//
//    }
//    catch (const exception& e) {
//        cerr << "Error during simulation: " << e.what() << endl;
//    }
//
//    return 0;
//}








