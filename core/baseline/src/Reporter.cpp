//#include "Reporter.h"
//#include <iostream>
//#include <fstream>
//#include <iomanip> // For std::fixed, std::setprecision, and std::setw
//
//using namespace BaseLine;
//using namespace std;
//
//void Reporter::report(const std::string& filename, const std::vector<Coords3D>& positions, const std::vector<Coords3D>& velocities, const std::vector<Coords3D>& forces, int step, bool clearFile) {
//    // Determine file open mode based on whether clearing the file or appending
//    ios_base::openmode fileMode = clearFile ? (ios::out) : (ios::out | ios::app);
//
//    // Open the file with the selected mode
//    ofstream outputFile(filename, fileMode);
//
//    if (!outputFile.is_open()) {
//        cerr << "Unable to open the file: " << filename << endl;
//        return;
//    }
//
//    // Write headers only if clearing the file
//    if (clearFile) {
//        outputFile << "Step\tID\tPosition\t\t\t\t\tVelocity\t\t\t\t\tForce\n";
//    }
//
//    // Set fixed number of decimal places for floating-point output
//    outputFile << fixed << setprecision(6);
//
//    int numAtoms = min(10, static_cast<int>(positions.size()));
//    for (int i = 0; i < numAtoms; ++i) {
//        outputFile << step << "\t" << i + 1
//            << "\t(" << setw(10) << positions[i][0] << ", " << setw(10) << positions[i][1] << ", " << setw(10) << positions[i][2] << ")\t"
//            << "(" << setw(10) << velocities[i][0] << ", " << setw(10) << velocities[i][1] << ", " << setw(10) << velocities[i][2] << ")\t"
//            << "(" << setw(10) << forces[i][0] << ", " << setw(10) << forces[i][1] << ", " << setw(10) << forces[i][2] << ")\n";
//    }
//
//    outputFile.close();
//}


#include "stdafx.h"
#include "Reporter.h"
#include <iomanip> // For std::fixed, std::setprecision, and std::setw

using namespace std;
using namespace BaseLine;


//std::ofstream Reporter::outputFile;
Reporter::Reporter() {
    // Constructor implementation (if necessary, otherwise leave empty)
}
void Reporter::TestPVFReport(const std::string& filename, const std::vector<Coords3D>& positions,
    const std::vector<Coords3D>& velocities, const std::vector<Coords3D>& forces, int step, const vector<PTorsionParams>& torsionParams, const vector<HBondParams>& bondParams, const vector<HAngleParams>& angleParams) {


    // Open the file with the appropriate mode
    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(filename, fileMode);

    //ofstream outputFile(filename, ios::out | ((step == 0) ? ios::trunc : ios::app));

    //if (!outputFile.is_open()) {
    //    throw runtime_error("Unable to open file: " + filename);
    //}

    //outputFile << fixed << setprecision(3);



    // Check if the file is successfully opened
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename << endl;
        return;
    }

    //if (step == 0) {
    //    outputFile << "Step\tID\tPosition (nm)\t\tVelocity (nm/ps)\t\tForce (kJ/mol*nm)\n";
    //}

    outputFile << "Step" << step<<"\n";
    outputFile << "ID\tPosition(nm)\t\t\tVelocity(nm / ps)\t\tForce(kJ / mol * nm)\n";

    outputFile << fixed << setprecision(3);

    //int numAtoms = min(10, static_cast<int>(positions.size()));
    //for (int i = 0; i < numAtoms; ++i) {
    //    outputFile
    //        << i + 1
    //        << "\t(" << setw(6) << positions[i][0] << ", " << setw(6) << positions[i][1] << ", " << setw(6) << positions[i][2] << ")\t"
    //        << "(" << setw(6) << velocities[i][0] << ", " << setw(6) << velocities[i][1] << ", " << setw(6) << velocities[i][2] << ")\t"
    //        << "(" << setw(7) << forces[i][0] << ", " << setw(7) << forces[i][1] << ", " << setw(7) << forces[i][2] << ")\n";


    //}

    //vector<int> atoms;

    //if (!torsionParams.empty()) {

    //    //For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0] and only print out the first torsion atom ids
    //    atoms = { torsionParams[0].p1, torsionParams[0].p2, torsionParams[0].p3, torsionParams[0].p4 };
    //    //cout << atoms[0]<<"  " << atoms[1] << "  " << atoms[2] << "  " << atoms[3] << "  ";
    //}
    //if (!bondParams.empty()) {

    //    //For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0] and only print out the first torsion atom ids
    //    atoms = { bondParams[0].p1, bondParams[0].p2 };
    //    //cout << atoms[0]<<"  " << atoms[1] << "  " << atoms[2] << "  " << atoms[3] << "  ";
    //}
    //if (!angleParams.empty()) {

    //    //For debugging purposes and PTF diagnosis, I only look at the first torsion so torsionParams[0] and only print out the first torsion atom ids
    //    atoms = { angleParams[0].p1, angleParams[0].p2, angleParams[0].p3};
    //    //cout << atoms[0]<<"  " << atoms[1] << "  " << atoms[2] << "  " << atoms[3] << "  ";
    //}


    for (int i = 0; i < positions.size(); ++i) {  // Loop over the 4 atoms in the first torsion
        //int atomIndex = atoms[i];  // Get the actual atom index
        outputFile
            << i + 1
            //for double precision checking
            //<< "\t(" << setw(6) << positions[i][0] << ", " << setw(6) << positions[i][1] << ", " << setw(6) << positions[i][2] << ")\t"
            //<< "(" << setw(6) << velocities[i][0] << ", " << setw(6) << velocities[i][1] << ", " << setw(6) << velocities[i][2] << ")\t"
            //<< setprecision(20) << "(" << setw(24) << forces[i][0] << ", " << setw(24) << forces[i][1] << ", " << setw(24) << forces[i][2] << ")\n" << setprecision(3);
            //for a cleaner output
            << "\t(" << setprecision(10) << setw(10) << positions[i][0] << ", " << setw(10) << positions[i][1] << ", " << setw(10) << positions[i][2] << ")\t"
            << "(" << setw(10) << velocities[i][0] << ", " << setw(10) << velocities[i][1] << ", " << setw(10) << velocities[i][2] << ")\t"
             << "(" << setw(10) << forces[i][0] << ", " << setw(10) << forces[i][1] << ", " << setw(10) << forces[i][2] << ")\n" ;
    }



    outputFile << "\n";
    outputFile.close();

}





//     static void TotalEnergy(const std::string& baseFilename, const std::vector<double>& kineticEnergies, double totalPotentialEnergy, int step);



//    void TotalEnergyReport(const std::string& baseFilename, const std::vector<double>& kineticEnergies, double totalPotentialEnergy, int step);


void Reporter::TotalEnergyReport(const string& baseFilename, double& totalKEnergy, double& totalPEnergy, double& totalEnergy, int step) {
    std::string filename = "TotalEnergy" + baseFilename;
    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(filename, fileMode);

    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename << endl;
        return;
    }

    if (step == 0) {
        outputFile << fixed << setprecision(10);
        outputFile << "Step\tTotal Kinetic Energy (kJ/mol)\tTotal Potential Energy (kJ/mol)\tTotal Energy (kJ/mol)\n";
    }

    outputFile << setw(5) << step
        << setw(25) << totalKEnergy
        << setw(30) << totalPEnergy
        << setw(25) << totalEnergy << "\n";

    outputFile.close();
}