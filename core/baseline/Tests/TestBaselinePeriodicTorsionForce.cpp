#include "stdafx.h"
#include "Engine.h"
#include "Coords3D.h"

using namespace std;
using namespace BaseLine;

int main() {

    // Assuming these vectors need to be passed to the RunSimulation method.
    vector<Coords3D> totalForces; // Initialized empty, assuming Engine handles their setup.
    vector<Coords3D> velocities; // Initialized empty, similarly.


    // first torsion info: 
    // <Torsion k=".6508444444444444" p1="0" p2="4" p3="6" p4="7" periodicity="3" phase="0"/>

    vector<Coords3D> atomPositions = {
        {2.5160000324249268, 1.4160000085830688, 1.9440000057220459},
        {2.5090000629425049, 1.3919999599456787, 1.7979999780654907},
        {2.4419999122619629, 1.4969999790191650, 1.7100000381469727},
        {2.4590001106262207, 1.4630000591278076, 1.6080000400543213}
    };

    vector<PTorsionParams> torsionParams = { {0, 1, 2, 3, 0.6508444444444444, 0.0, 3} };//order p1, p2, p3, p4, k double, phase double, periodicity int
    vector<HAngleParams> angleParams{};//empty as in here we only want to test the accuracy of software result for 1 periodic torsion force
    vector<double> masses = { 14.01, 12.01, 12.01, 1.008 };

    Engine engine("","",atomPositions, masses, torsionParams, angleParams);//engine construction for test run
    string outputFilename = "output_test_PTF_1Torsion.txt"; // Output file for reporting

    double StepSize = 0.001; // Example step size
    int TotalSteps = 2; // Example total number of steps for the simulation

    // Run the simulation with specified parameters
    try {
        engine.RunSimulation(outputFilename, StepSize, TotalSteps);// , systemFilename, stateFilename, totalForces, velocities);

        cout << "Simulation completed successfully." << endl;
    }
    catch (const exception& e) {
        cerr << "Error during simulation: " << e.what() << endl;
    }

    return 0;
}