//#include "stdafx.h"
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "Reporter.h"
#include <iomanip> // For fixed, setprecision, and setw
#include "Coords3D.h" 
#include <chrono>
#include <cctype>  // For isupper and islower



using namespace std;
using namespace BaseLine;


//ofstream Reporter::outputFile;
Reporter::Reporter() {
    // Constructor implementation (if necessary, otherwise leave empty)
}
void Reporter::TestPVFReport(const string& filename, const vector<Coords3D>& positions,
    const vector<Coords3D>& velocities, const vector<Coords3D>& forces, int step, const vector<PTorsionParams>& torsionParams, const vector<BondParams>& bondParams, const vector<AngleParams>& angleParams) {
    string filename_modified = "PVFReport" + filename;


    // Open the file with the appropriate mode
    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(filename_modified, fileMode);

    //ofstream outputFile(filename_modified, ios::out | ((step == 0) ? ios::trunc : ios::app));

    //if (!outputFile.is_open()) {
    //    throw runtime_error("Unable to open file: " + filename_modified);
    //}

    //outputFile << fixed << setprecision(3);



    // Check if the file is successfully opened
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename_modified << endl;
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





//     static void TotalEnergy(const string& baseFilename, const vector<double>& kineticEnergies, double totalPotentialEnergy, int step);



//    void TotalEnergyReport(const string& baseFilename, const vector<double>& kineticEnergies, double totalPotentialEnergy, int step);


void Reporter::TotalEnergyReport(const string& baseFilename, double& totalKEnergy, double& totalPEnergy, double& totalEnergy, int step) {
    string filename = "TotalEnergy" + baseFilename;
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


// Helper function to get the current date in YYYY-MM-DD format
string Reporter::getCurrentDate() {
    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);
    stringstream ss;
    ss << put_time(localtime(&in_time_t), "%Y-%m-%d");
    return ss.str();
}

// Extracts the element from the atom name
string Reporter::extractElement(const string& atomName) {
    if (atomName.empty()) return "";
    string element;
    element += atomName[0]; // First character is always included and is uppercase
    if (atomName.length() > 1 && islower(atomName[1])) {
        // Include the second character only if it is lowercase
        element += atomName[1];
    }
    return element;
}

PDBAtom Reporter::parsePDBLine(const string& line) {
    PDBAtom atom;
    stringstream ss(line);
    ss >> atom.recordType >> atom.atomNum >> atom.atomName >> atom.resName >> atom.resSeq
        >> atom.coords[0] >> atom.coords[1] >> atom.coords[2] >> atom.occupancy >> atom.tempFactor;
    atom.chainID = 'A';//this is a simplification for now
    atom.element = extractElement(atom.atomName);
    atom.recordType = (standardAminoAcids.count(atom.resName) > 0) ? "ATOM  " : "HETATM";
    return atom;
}



void Reporter::outputTemplateExtractor(const string& inputFilename, vector<PDBAtom>& outputTemplate) {
    ifstream inputFile(inputFilename);
    if (!inputFile) {
        cerr << "Failed to open input file." << endl;
        return;
    }
    string line;
    while (getline(inputFile, line)) {
        if (line.substr(0, 4) == "ATOM") {
            PDBAtom atom = parsePDBLine(line);
            outputTemplate.push_back(atom);
        }
        else if (line.substr(0, 6) == "CRYST1") {
            _PBCLine = line;
        }

    }
    inputFile.close();
}


// Helper function to format a PDB line
//string Reporter::formatPDBLine(const PDBAtom& atom) {
//    stringstream ss;
//    //ss << setw(6) << left << atom.recordType
//    //    << setw(5) << right << atom.atomNum << " "
//    //    << setw(4) << left << atom.atomName
//    //    << setw(3) << right << atom.resName << " "
//    //    << atom.chainID
//    //    << setw(4) << right << atom.resSeq << "    "
//    //    << fixed << setprecision(3)
//    //    << setw(8) << right << atom.coords[0]
//    //    << setw(8) << right << atom.coords[1]
//    //    << setw(8) << right << atom.coords[2]
//    //    << setw(6) << atom.occupancy
//    //    << setw(6) << atom.tempFactor << "          "
//    //    << setw(2) << left << atom.element;
//    //return ss.str();
//
//    ss << setw(6) << left << atom.recordType << " "
//        << setw(5) << right << atom.atomNum << "  "
//        << setw(5) << left << atom.atomName << " "
//        << setw(3) << right << atom.resName << " "
//        << atom.chainID << " "
//        << setw(4) << right << atom.resSeq << "    "
//        << fixed << setprecision(3)
//        << setw(8) << left << atom.coords[0]
//        << setw(8) << left << atom.coords[1]
//        << setw(8) << left << atom.coords[2]<<"  "
//        << fixed << setprecision(2)
//        << setw(6) << left << atom.occupancy << "  "
//        << setw(6) << left << atom.tempFactor << "          "
//        << setw(2) << left << atom.element;
//    return ss.str();
//
//}


// Helper function to format a PDB line
string Reporter::formatPDBLine(const PDBAtom& atom) {
    stringstream ss;

    // Format string according to PDB column alignment
    ss << left << setw(6) << atom.recordType
        << right << setw(5) << atom.atomNum << " "
        << left << setw(4) << atom.atomName
        << right << setw(1) << " " // AltLoc placeholder
        << setw(3) << atom.resName << " "
        << setw(1) << atom.chainID
        << setw(4) << right << atom.resSeq
        << setw(1) << " " // iCode placeholder
        << "   " // three spaces before coordinates start
        << fixed << setprecision(3)
        << right << setw(8) << atom.coords[0]
        << right << setw(8) << atom.coords[1]
        << right << setw(8) << atom.coords[2]
        << fixed << setprecision(2)
        << setw(6) << atom.occupancy
        << setw(6) << atom.tempFactor
        << "          " // 10 spaces before element
        << left << setw(2) << atom.element;

    return ss.str();
}



// Main function to read, process, and write PDB file
//part one is triggered when step=0 and the output template (outputTemplate) is going to be constructed
void Reporter::pdbOutputGeneratorPart1(const string& inputFilename, const string& outputFilename, vector<PDBAtom>& outputTemplate) {

    //ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(outputFilename);


    outputTemplateExtractor(inputFilename, outputTemplate);//the outputTemplate only needs to be generated once
    outputFile << "REMARK   Generated by NexaBind 1.0, " << getCurrentDate() << endl;
    outputFile << _PBCLine << endl;

}


// Main function to read, process, and write PDB file
void Reporter::pdbOutputGeneratorPart2(const string& outputFilename, vector<PDBAtom>& outputTemplate, const vector<Coords3D>& positions, int ModelNum) {

    ios_base::openmode fileMode = (ios::out | ios::app);
    ofstream outputFile(outputFilename, fileMode);


    // Check if the file is successfully opened
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << outputFilename << endl;
        return;
    }

    // Write the REMARK line with the current date
    // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)
    cout << ModelNum << "   ";
    
    outputFile << "MODEL        " << ModelNum  << endl;//in pdb step starts from 1


    for (int i = 0; i < outputTemplate.size(); ++i) {
        const auto& atom = outputTemplate[i];
        PDBAtom updatedAtom = atom;
        // bug fixed: positions are in nm and need to be multiplied to 10 to be in A unit as the output.pdb should be in A
        updatedAtom.coords = positions[atom.atomNum - 1]*10.0; // Adjust coordinates
        outputFile << formatPDBLine(updatedAtom) << endl;

        // Check if the current atom is the last one in the vector
        if (i == outputTemplate.size() - 1) {
            //outputFile << setw(6) << left << "TER" << " " << setw(5) << right << atom.atomNum + 1
            //    << setw(8) << ""
            //    << setw(3) << left << atom.resName << " "
            //    << atom.chainID << " "
            //    << setw(4) << right << atom.resSeq << endl;


            outputFile << left << setw(6) << "TER"
                << right << setw(5) << atom.atomNum + 1 << " "
                << left << setw(4) << ""
                << right << setw(1) << " " // AltLoc placeholder
                << setw(3) << atom.resName << " "
                << setw(1) << atom.chainID
                << setw(4) << right << atom.resSeq << endl;

        }

    }

    outputFile << "ENDMDL" << endl;
    outputFile.close();
}

//// Main function to read, process, and write PDB file
//void Reporter::pdbOutputGenerator(const string& inputFilename, const string& outputFilename, vector<PDBAtom>& outputTemplate, const vector<Coords3D>& positions, int step) {
//
//    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
//    ofstream outputFile(outputFilename, fileMode);
//
//    if (step == 0) {
//        outputTemplateExtractor(inputFilename, outputTemplate);//the outputTemplate only needs to be generated once
//        outputFile << "REMARK   Generated by NexaBind 1.0, " << getCurrentDate() << endl;
//        outputFile << _PBCLine << endl;
//    }
//
//    // Check if the file is successfully opened
//    if (!outputFile.is_open()) {
//        cerr << "Unable to open the file: " << outputFilename << endl;
//        return;
//    }
//
//    // Write the REMARK line with the current date
//    // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)
//    outputFile << "MODEL        " << step + 1 << endl;//in pdb step starts from 1
//
//
//    for (size_t i = 0; i < outputTemplate.size(); ++i) {
//        const auto& atom = outputTemplate[i];
//        PDBAtom updatedAtom = atom;
//        updatedAtom.coords = positions[atom.atomNum - 1]; // Adjust coordinates
//        outputFile << formatPDBLine(updatedAtom) << endl;
//
//        // Check if the current atom is the last one in the vector
//        if (i == outputTemplate.size() - 1) {
//            outputFile << setw(6) << left << "TER" << " " << setw(5) << right << atom.atomNum + 1
//                << setw(8) << ""
//                << setw(3) << left << atom.resName << " "
//                << atom.chainID << " "
//                << setw(4) << right << atom.resSeq << endl;
//        }
//    }
//
//    outputFile << "ENDMDL" << endl;
//    outputFile.close();
//}

