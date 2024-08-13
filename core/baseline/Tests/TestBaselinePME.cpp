#pragma once

//#include <string>
//#include <sstream>
//#include <fstream>
//#include <utility>
#include <iostream>
#include <vector>
#include "Engine.h"
#include "Coords3D.h"

using namespace std;
using namespace BaseLine;

int main() {
    vector<Coords3D> atomPositions;
    vector<double> masses;

    //vector<PTorsionParams> torsionParams = {};//order p1, p2, p3, p4, k double, phase double, periodicity int
    //vector<HBondParams> bondParams{};//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant
    //vector<HAngleParams> angleParams{};//empty as in here we only want to test the accuracy of software result for 1 periodic torsion force

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


    vector<HAngleParams> angleParams = {// this is added here purely for creating exclusions
    {0, 1, 2, 1.911135530933791, 292.88},     // p1, p2, p3, angle, k
    {0, 1, 3, 1.911135530933791, 292.88},     // p1, p2, p3, angle, k
    {0, 1, 4, 1.911135530933791, 418.40000000000003},     // p1, p2, p3, angle, k
    {1, 4, 5, 2.101376419401173, 669.44},     // p1, p2, p3, angle, k
    {1, 4, 6, 2.035053907825388, 585.76},     // p1, p2, p3, angle, k
    {2, 1, 3, 1.911135530933791, 292.88},     // p1, p2, p3, angle, k
    {2, 1, 4, 1.911135530933791, 418.40000000000003},     // p1, p2, p3, angle, k
    {3, 1, 4, 1.911135530933791, 418.40000000000003},     // p1, p2, p3, angle, k
    {4, 6, 7, 2.0943951023931953, 418.40000000000003},    // p1, p2, p3, angle, k
    {4, 6, 8, 2.1275563581810877, 418.40000000000003},    // p1, p2, p3, angle, k
    {5, 4, 6, 2.1450096507010312, 669.44},    // p1, p2, p3, angle, k
    {6, 8, 9, 1.911135530933791, 418.40000000000003},     // p1, p2, p3, angle, k
    {6, 8, 10, 1.9146261894377796, 669.44},   // p1, p2, p3, angle, k
    {6, 8, 14, 1.9216075064457567, 527.184},  // p1, p2, p3, angle, k
    {7, 6, 8, 2.0601866490541068, 418.40000000000003},    // p1, p2, p3, angle, k
    {8, 10, 11, 1.911135530933791, 418.40000000000003},   // p1, p2, p3, angle, k
    {8, 10, 12, 1.911135530933791, 418.40000000000003},   // p1, p2, p3, angle, k
    {8, 10, 13, 1.911135530933791, 418.40000000000003},   // p1, p2, p3, angle, k
    {8, 14, 15, 2.101376419401173, 669.44},   // p1, p2, p3, angle, k
    {8, 14, 16, 2.035053907825388, 585.76},   // p1, p2, p3, angle, k
    {9, 8, 10, 1.911135530933791, 418.40000000000003},    // p1, p2, p3, angle, k
    {9, 8, 14, 1.911135530933791, 418.40000000000003},    // p1, p2, p3, angle, k
    {10, 8, 14, 1.9390607989657, 527.184},    // p1, p2, p3, angle, k
    {11, 10, 12, 1.911135530933791, 292.88},  // p1, p2, p3, angle, k
    {11, 10, 13, 1.911135530933791, 292.88},  // p1, p2, p3, angle, k
    {12, 10, 13, 1.911135530933791, 292.88},  // p1, p2, p3, angle, k
    {14, 16, 17, 2.0943951023931953, 418.40000000000003}, // p1, p2, p3, angle, k
    {14, 16, 18, 2.1275563581810877, 418.40000000000003}, // p1, p2, p3, angle, k
    {15, 14, 16, 2.1450096507010312, 669.44}, // p1, p2, p3, angle, k
    {16, 18, 19, 1.911135530933791, 418.40000000000003},  // p1, p2, p3, angle, k
    {16, 18, 20, 1.911135530933791, 418.40000000000003},  // p1, p2, p3, angle, k
    {16, 18, 21, 1.911135530933791, 418.40000000000003},  // p1, p2, p3, angle, k
    {17, 16, 18, 2.0601866490541068, 418.40000000000003}, // p1, p2, p3, angle, k
    {19, 18, 20, 1.911135530933791, 292.88},  // p1, p2, p3, angle, k
    {19, 18, 21, 1.911135530933791, 292.88},  // p1, p2, p3, angle, k
    {20, 18, 21, 1.911135530933791, 292.88}   // p1, p2, p3, angle, k
    };




    vector<PTorsionParams> torsionParams = {// this is added here purely for creating exclusions
    {0, 1, 4, 5, 3.3472, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {0, 1, 4, 5, 0.33472, 3.141592653589793, 3},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {1, 4, 6, 7, 10.46, 3.141592653589793, 2},       // p1, p2, p3, p4, k double, phase double, periodicity int
    {1, 4, 6, 8, 10.46, 3.141592653589793, 2},       // p1, p2, p3, p4, k double, phase double, periodicity int
    {2, 1, 4, 5, 3.3472, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {2, 1, 4, 5, 0.33472, 3.141592653589793, 3},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {3, 1, 4, 5, 3.3472, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {3, 1, 4, 5, 0.33472, 3.141592653589793, 3},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 6, 8, 10, 3.3472, 0.0, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 6, 8, 10, 7.5312, 0.0, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 6, 8, 10, 8.368, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 6, 8, 14, 1.75728, 0.0, 3},   // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 6, 8, 14, 1.12968, 0.0, 2},   // p1, p2, p3, p4, k double, phase double, periodicity int
    {5, 4, 6, 7, 10.46, 3.141592653589793, 2},       // p1, p2, p3, p4, k double, phase double, periodicity int
    {5, 4, 6, 7, 8.368, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {5, 4, 6, 8, 10.46, 3.141592653589793, 2},       // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 10, 11, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 10, 12, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 10, 13, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 14, 16, 2.3012, 3.141592653589793, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 14, 16, 6.610720000000001, 3.141592653589793, 2}, // p1, p2, p3, p4, k double, phase double, periodicity int
    {6, 8, 14, 16, 1.8828, 3.141592653589793, 1},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {8, 14, 16, 17, 10.46, 3.141592653589793, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {8, 14, 16, 18, 10.46, 3.141592653589793, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {9, 8, 10, 11, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {9, 8, 10, 12, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {9, 8, 10, 13, 0.6508444444444444, 0.0, 3},      // p1, p2, p3, p4, k double, phase double, periodicity int
    {9, 8, 14, 15, 3.3472, 0.0, 1},     // p1, p2, p3, p4, k double, phase double, periodicity int
    {9, 8, 14, 15, 0.33472, 3.141592653589793, 3},   // p1, p2, p3, p4, k double, phase double, periodicity int
    {10, 8, 14, 16, 1.6736, 0.0, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {10, 8, 14, 16, 0.8368, 0.0, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {10, 8, 14, 16, 0.8368, 0.0, 1},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {11, 10, 8, 14, 0.6508444444444444, 0.0, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {12, 10, 8, 14, 0.6508444444444444, 0.0, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {13, 10, 8, 14, 0.6508444444444444, 0.0, 3},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {15, 14, 16, 17, 10.46, 3.141592653589793, 2},   // p1, p2, p3, p4, k double, phase double, periodicity int
    {15, 14, 16, 17, 8.368, 0.0, 1},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {15, 14, 16, 18, 10.46, 3.141592653589793, 2},   // p1, p2, p3, p4, k double, phase double, periodicity int
    {1, 6, 4, 5, 43.932, 3.141592653589793, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {4, 8, 6, 7, 4.6024, 3.141592653589793, 2},    // p1, p2, p3, p4, k double, phase double, periodicity int
    {8, 16, 14, 15, 43.932, 3.141592653589793, 2},  // p1, p2, p3, p4, k double, phase double, periodicity int
    {14, 18, 16, 17, 4.6024, 3.141592653589793, 2}  // p1, p2, p3, p4, k double, phase double, periodicity int
    };





    // Initialize nonbonded parameters
    NonbondedParams nbParams;
    nbParams.alpha = 0.3; // Example value for Ewald summation alpha parameter
    nbParams.cutoff = 1.0; // Nonbonded cutoff distance in nanometers
    nbParams.ewaldTolerance = 0.0005; // Tolerance for Ewald summation
    nbParams.rfDielectric = 78.5; // Reaction field dielectric constant
    nbParams.switchingDistance = 0.9; // Switching distance for force switching
    nbParams.dispersionCorrection = 1;
    nbParams.exceptionsUsePeriodic = 1;
    nbParams.forceGroup = 0;
    nbParams.includeDirectSpace = 1;
    nbParams.ljAlpha = 2;
    nbParams.ljnx = 10;
    nbParams.ljny = 10;
    nbParams.ljnz = 10;
    nbParams.method = 2; // PME method identifier, assuming 2 stands for PME
    nbParams.nx = 32;
    nbParams.ny = 32;
    nbParams.nz = 32;
    nbParams.recipForceGroup = 1;
    nbParams.useSwitchingFunction = 0;
    nbParams.version = 1;

    string pme_force_test = "Dialanine";// another option is NACl
    if (pme_force_test == "Dialanine") {
        const int numParticles = 22;
        //double masses[] = { 1.008, 12.01, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 12.01, 1.008, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 1.008, 1.008 };  // Specific masses for the particles

        masses = { 1.008, 12.01, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 12.01, 1.008, 1.008, 1.008, 12.01, 16, 14.01, 1.008, 12.01, 1.008, 1.008, 1.008 };  // Specific masses for the particles


        //// Set the relative dielectric constant of water
        //double relativeDielectricConstant = 78.54;
        //nonbondedForce->setReactionFieldDielectric(relativeDielectricConstant);



        // Add particles with properties
        nbParams.particles.push_back({ 0.1123, 0.2649532787749369, 0.06568879999999999 });// how reads this: charge, sigma, epsilon
        nbParams.particles.push_back({ -0.3662, 0.3399669508423535, 0.4577296 });
        nbParams.particles.push_back({ 0.1123, 0.2649532787749369, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.1123, 0.2649532787749369, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.5972, 0.3399669508423535, 0.359824 });
        nbParams.particles.push_back({ -0.5679, 0.2959921901149463, 0.87864 });
        nbParams.particles.push_back({ -0.4157, 0.3249998523775958, 0.7112800000000001 });
        nbParams.particles.push_back({ 0.2719, 0.10690784617684071, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.0337, 0.3399669508423535, 0.4577296 });
        nbParams.particles.push_back({ 0.0823, 0.2471353044121301, 0.06568879999999999 });
        nbParams.particles.push_back({ -0.1825, 0.3399669508423535, 0.4577296 });
        nbParams.particles.push_back({ 0.0603, 0.2649532787749369, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.0603, 0.2649532787749369, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.0603, 0.2649532787749369, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.5973, 0.3399669508423535, 0.359824 });
        nbParams.particles.push_back({ -0.5679, 0.2959921901149463, 0.87864 });
        nbParams.particles.push_back({ -0.4157, 0.3249998523775958, 0.7112800000000001 });
        nbParams.particles.push_back({ 0.2719, 0.10690784617684071, 0.06568879999999999 });
        nbParams.particles.push_back({ -0.149, 0.3399669508423535, 0.4577296 });
        nbParams.particles.push_back({ 0.0976, 0.2471353044121301, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.0976, 0.2471353044121301, 0.06568879999999999 });
        nbParams.particles.push_back({ 0.0976, 0.2471353044121301, 0.06568879999999999 });






        //I ignore the exceptions for now
        //// Add exceptions, if any  // i don't include them for now since NexaBind doesn't have exceptions
        //nonbondedForce->addException(0, 1, 0, 1, 0);
        //nonbondedForce->addException(0, 2, 0, 1, 0);
        //nonbondedForce->addException(1, 2, 0, 1, 0);
        //nonbondedForce->addException(0, 3, 0, 1, 0); ....
        // ...


        atomPositions = {
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
    }
    else if (pme_force_test == "NACl") {


        // Setup for Na+ and Cl- ions
        atomPositions = {

            {1.0, 1.0, 1.0}, // Position for Na+
            {2.0, 2.0, 2.0}  // Position for Cl-
        };

        // Masses for Na+ and Cl-, in atomic mass units
        masses = { 22.989769, 35.453 };

        // Add nonbonded particles, assuming the class has a method to handle this
        nbParams.particles.push_back({ 1.0, 0.235, 0.352 * 4.184 });  // Na+
        nbParams.particles.push_back({ -1.0, 0.4401, 0.100 * 4.184 });  // Cl-

    }



    PeriodicBoundaryCondition::BoxInfo boxInfo;
    boxInfo.boxSize = { 3.0, 3.0, 3.0 }; // Example manual setting // note that 
    // Initialize the engine with nonbonded parameters
    Engine engine("", "", atomPositions, masses, torsionParams, bondParams, angleParams, nbParams, boxInfo);
    //Engine engine("", "", atomPositions, masses, {}, {}, {}, nbParams, boxInfo);

    string outputFilename = "output_PME_test.txt"; // Output file for reporting

    double StepSize = 0.001; // Example step size in picoseconds
    int TotalSteps = 10; // Example total number of steps for the simulation
    int interval = 1; // The interval at which the reporters will save data

    // Run the simulation with specified parameters
    try {
        engine.RunSimulation("", outputFilename, StepSize, TotalSteps, interval);
        cout << "PME test simulation completed successfully." << endl;
    }
    catch (const exception& e) {
        cerr << "Error during simulation: " << e.what() << endl;
    }

    return 0;
}
