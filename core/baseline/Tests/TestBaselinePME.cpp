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

    vector<PTorsionParams> torsionParams = {};//order p1, p2, p3, p4, k double, phase double, periodicity int
    vector<HBondParams> bondParams{};//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant
    vector<HAngleParams> angleParams{};//empty as in here we only want to test the accuracy of software result for 1 periodic torsion force


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
