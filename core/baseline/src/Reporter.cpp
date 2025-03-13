//#include "stdafx.h"
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include "Reporter.h"
#include <iomanip> 
#include "Coords3D.h" 
#include <chrono>
#include <cctype> 



using namespace std;
using namespace BaseLine;


//ofstream Reporter::outputFile;
Reporter::Reporter() {
    // Constructor implementation (if necessary, otherwise leave empty)
}
void Reporter::TestPVFReport(const string& filename, const vector<Coords3D>& positions,
    const vector<Coords3D>& velocities, const vector<Coords3D>& forces, int step, const vector<PTorsionParams>& torsionParams, const vector<BondParams>& bondParams, const vector<AngleParams>& angleParams) {
    string filename_modified = "PVFReport" + filename;


    // To open the file with the appropriate mode
    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(filename_modified, fileMode);


    // To check if the file is successfully opened
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename_modified << endl;
        return;
    }

    outputFile << "Step" << step<<"\n";
    outputFile << "ID\tPosition(nm)\t\t\tVelocity(nm / ps)\t\tForce(kJ / mol * nm)\n";

    outputFile << fixed << setprecision(3);



    for (int i = 0; i < positions.size(); ++i) {  // Loop over the 4 atoms in the first torsion
        outputFile
            << i + 1
            << "\t(" << setprecision(10) << setw(10) << positions[i][0] << ", " << setw(10) << positions[i][1] << ", " << setw(10) << positions[i][2] << ")\t"
            << "(" << setw(10) << velocities[i][0] << ", " << setw(10) << velocities[i][1] << ", " << setw(10) << velocities[i][2] << ")\t"
             << "(" << setw(10) << forces[i][0] << ", " << setw(10) << forces[i][1] << ", " << setw(10) << forces[i][2] << ")\n" ;
    }



    outputFile << "\n";
    outputFile.close();

}

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
string Reporter::formatPDBLine(const PDBAtom& atom) {
    stringstream ss;

    // To format string according to PDB column alignment
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
    outputFile << "REMARK   Generated by Nexaus 1.0.0, " << getCurrentDate() << endl;
    outputFile << _PBCLine << endl;

}


// Main function to read, process, and write PDB file
void Reporter::pdbOutputGeneratorPart2(const string& outputFilename, vector<PDBAtom>& outputTemplate, const vector<Coords3D>& positions, int ModelNum) {

    ios_base::openmode fileMode = (ios::out | ios::app);
    ofstream outputFile(outputFilename, fileMode);


    // To check if the file is successfully opened
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << outputFilename << endl;
        return;
    }

    // To write the REMARK line with the current date
    // Model number is different than step number. it's step divided by interval (the interval at which the reporters will save data)
    cout << ModelNum << "   ";
    
    outputFile << "MODEL        " << ModelNum  << endl;//in pdb step starts from 1


    for (int i = 0; i < outputTemplate.size(); ++i) {
        const auto& atom = outputTemplate[i];
        PDBAtom updatedAtom = atom;
        updatedAtom.coords = positions[atom.atomNum - 1]*10.0; // Adjust coordinates
        outputFile << formatPDBLine(updatedAtom) << endl;

        // To check if the current atom is the last one in the vector
        if (i == outputTemplate.size() - 1) {
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

// Function to report thermostat and barostat-related properties
void Reporter::EnsembleReport(const string& baseFilename, double& totalKEnergy, double& calcTemp, double& volume, double& density, int& step) {
    // Open file in append mode if step > 0, otherwise create a new file

    string filename = "EnsembleReport" + baseFilename;
    ios_base::openmode fileMode = (step == 0) ? (ios::out) : (ios::out | ios::app);
    ofstream outputFile(filename, fileMode);

    // Check if the file opened successfully
    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename << endl;
        return;
    }

    // Write header if this is the first step
    if (step == 0) {
        outputFile << "#Step,Kinetic Energy (kJ/mole),Temperature (K),Box Volume (nm^3),Density (g/mL)\n";
    }

    // Append data for the current step
    outputFile << step << ",  "
        << fixed << setprecision(6) << totalKEnergy << ",  "
        << fixed << setprecision(6) << calcTemp << ",  "
        << fixed << setprecision(6) << volume << ",  "
        << fixed << setprecision(6) << density << "  \n";

    // Close file
    outputFile.close();
}