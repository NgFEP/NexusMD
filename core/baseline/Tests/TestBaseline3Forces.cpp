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


    // first torsion info: 
    // <Torsion k=".6508444444444444" p1="0" p2="4" p3="6" p4="7" periodicity="3" phase="0"/>

    //vector<Coords3D> atomPositions = {
    //    {2.5160000324249268, 1.4160000085830688, 1.9440000057220459},
    //    {2.5090000629425049, 1.3919999599456787, 1.7979999780654907},
    //    {2.4419999122619629, 1.4969999790191650, 1.7100000381469727},
    //    {2.4590001106262207, 1.4630000591278076, 1.6080000400543213}
    //};


    // first harmonic bond force info: 
// <Bond d=".1522" k="265265.6" p1="19" p2="4"/>

    //vector<Coords3D> atomPositions = {
    //    {2.6519999504089355, 1.3969999551773071, 1.7430000305175781},//		<Position p1=19/>
    //    {2.509000062942505, 1.3919999599456787, 1.7979999780654907}//		<Position p2=4/>
    //};


    // first harmonic angle force info: 
// <Angle a="1.911135530933791" k="418.40000000000003" p1="0" p2="4" p3="5"/>

    //vector<Coords3D> atomPositions = {
    //    {2.5160000324249268, 1.4160000085830688, 1.9440000057220459},//		<Position p=0/>
    //    {2.509000062942505, 1.3919999599456787, 1.7979999780654907},//		<Position p=4/>
    //    {2.4700000286102295, 1.2929999828338623, 1.7760000228881836}//		<Position p=5/>

    //};




    // combined atom positions with making sure repetitive atoms combined into one by assigning the accurate ids to force parameters (torsionParams, bondParams and angleParams)  
    vector<Coords3D> atomPositions = {
        {2.5160000324249268, 1.4160000085830688, 1.9440000057220459},//		<Position p1=0/>
        {2.5090000629425049, 1.3919999599456787, 1.7979999780654907},//		<Position p2=4/>
        {2.4419999122619629, 1.4969999790191650, 1.7100000381469727},//		<Position p3=6/>
        {2.4590001106262207, 1.4630000591278076, 1.6080000400543213},//		<Position p4=7/>
        {2.6519999504089355, 1.3969999551773071, 1.7430000305175781},//		<Position p5=19/>
        {2.4700000286102295, 1.2929999828338623, 1.7760000228881836}//		<Position p6=5/>
    };


    vector<PTorsionParams> torsionParams = { {0, 1, 2, 3, 0.6508444444444444, 0.0, 3} };//order p1, p2, p3, p4, k double, phase double, periodicity int
    vector<HBondParams> bondParams{ {4, 1, .1522, 265265.6} };//p1,p2,d,k //double d; // Ideal bond distance and double k; // Force constant
    vector<HAngleParams> angleParams{ {0, 1, 5, 1.911135530933791, 418.40000000000003} };//p1,p2,p3,a,k //double a; // Ideal angle in radians and double k; // Force constant
    
    
    //vector<double> masses = { 14.01, 12.01, 12.01, 1.008 };
    //vector<double> masses = { 12.01, 12.01 };
    //vector<double> masses = { 14.01, 12.01, 1.008 };

    // combined masses
    //vector<double> masses_ids = { 0,      1,    2,     3,     4,     5 };
    vector<double> masses = { 14.01, 12.01, 12.01, 1.008, 12.01, 1.008 };




    Engine engine("", "", atomPositions, masses, torsionParams, bondParams, angleParams);//engine construction for test run
    string outputFilename = "output_test_3Forces.txt"; // Output file for reporting

    double StepSize = 0.001; // Example step size 1 fs
    int TotalSteps = 20; // Example total number of steps for the simulation
    int interval = 1; // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)

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