#include "stdafx.h"
#include "Engine.h"
#include "Coords3D.h"

using namespace std;
using namespace BaseLine;

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
    vector<double> masses = { 12.01000000000000, 12.01000000000000 };

    Engine engine("","",atomPositions, masses, torsionParams, bondParams, angleParams);//engine construction for test run
    string outputFilename = "output_test_HBF_1Bond.txt"; // Output file for reporting

    double StepSize = 0.001; // Example step size
    int TotalSteps = 3; // Example total number of steps for the simulation

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