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
        {2.318000078201294, 1.773800015449524, 0.6651999950408936},
        {2.2149999141693115, 1.7527999877929688, 0.6362000107765198},
        {2.1559998989105225, 1.7297999858856201, 0.7251999974250793},
        {2.203000068664551, 1.6748000383377075, 0.5612000226974487},
        {2.1619999408721924, 1.8797999620437622, 0.5712000131607056},
        {2.046999931335449, 1.8788000345230103, 0.5271999835968018},
        {2.242000102996826, 1.9817999601364136, 0.5382000207901001},
        {2.3350000381469727, 1.9758000373840332, 0.5771999955177307},
        {2.2079999446868896, 2.11680006980896, 0.498199999332428},
        {2.1089999675750732, 2.1477999687194824, 0.5332000255584717},
        {2.312999963760376, 2.2047998905181885, 0.5662000179290771},
        {2.4149999618530273, 2.180799961090088, 0.5371999740600586},
        {2.2899999618530273, 2.307800054550171, 0.5401999950408936},
        {2.296999931335449, 2.1958000659942627, 0.6732000112533569},
        {2.180999994277954, 2.1277999877929688, 0.3492000102996826},
        {2.180999994277954, 2.236799955368042, 0.29120001196861267},
        {2.1600000858306885, 2.0127999782562256, 0.2842000126838684},
        {2.1419999599456787, 1.9248000383377075, 0.3301999866962433},
        {2.128000020980835, 2.0237998962402344, 0.14319999516010284},
        {2.062000036239624, 1.9448000192642212, 0.10620000213384628},
        {2.072999954223633, 2.115799903869629, 0.12120000272989273},
        {2.2170000076293945, 2.0208001136779785, 0.08020000159740448}
    };

    vector<PTorsionParams> torsionParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force
    //vector<HBondParams> bondParams{ {0, 1, .1522000000000000, 265265.6000000000} };//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant

    vector<HBondParams> bondParams{// this is added here purely for creating exclusions
        {4, 1, 0.1522, 265265.6},     // p1, p2, d, k
        {4, 5, 0.1229, 476976.0},     // p1, p2, d, k
        {1, 0, 0.109, 284512.0},      // p1, p2, d, k
        {1, 2, 0.109, 284512.0},      // p1, p2, d, k
        {1, 3, 0.109, 284512.0},      // p1, p2, d, k
        {4, 6, 0.1335, 410032.0},     // p1, p2, d, k
        {14, 8, 0.1522, 265265.6},    // p1, p2, d, k
        {14, 15, 0.1229, 476976.0},   // p1, p2, d, k
        {8, 10, 0.1526, 259408.0},    // p1, p2, d, k
        {8, 9, 0.109, 284512.0},      // p1, p2, d, k
        {8, 6, 0.1449, 282001.6},     // p1, p2, d, k
        {10, 11, 0.109, 284512.0},    // p1, p2, d, k
        {10, 12, 0.109, 284512.0},    // p1, p2, d, k
        {10, 13, 0.109, 284512.0},    // p1, p2, d, k
        {7, 6, 0.101, 363171.2},      // p1, p2, d, k
        {14, 16, 0.1335, 410032.0},   // p1, p2, d, k
        {18, 19, 0.109, 284512.0},    // p1, p2, d, k
        {18, 20, 0.109, 284512.0},    // p1, p2, d, k
        {18, 21, 0.109, 284512.0},    // p1, p2, d, k
        {18, 16, 0.1449, 282001.6},   // p1, p2, d, k
        {17, 16, 0.101, 363171.2}     // p1, p2, d, k
    };

    vector<HAngleParams> angleParams = {};//empty as in here we only want to test the accuracy of software result for 1 harmonic bond force
    NonbondedParams nonbondedParams;// in this case nonbondedParams.particles = {} is automatically equal to zero

    vector<double> masses = { 1.008, 12.01, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 12.01, 1.008, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 1.008, 1.008 };  // Specific masses for the particles


    string inputFilename = "";
    string outputFilename = "output_test_HBF_1Bond.txt"; // Output file for reporting
    string systemFilename = "";
    string stateFilename = "";

    double stepSize = 0.001; // Example step size
    int totalSteps = 10; // Example total number of steps for the simulation
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








