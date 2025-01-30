#include "PDBResidueParser.h"

using namespace std;

// Constructor to initialize with file name

PDBResidueParser::PDBResidueParser() {};
PDBResidueParser::~PDBResidueParser() {};


// Utility function to trim whitespace
void PDBResidueParser::trim(std::string& str) {
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char c) { return !std::isspace(c); }));
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), str.end());
}

// Utility function to infer element from atom name
std::string PDBResidueParser::inferElementFromAtomName(const std::string& atomName) {
    if (atomName.empty()) return "";

    // To remove numerical characters from the atom name
    std::string filteredAtomName;
    for (char c : atomName) {
        if (!std::isdigit(c)) {
            filteredAtomName += c;
        }
    }

    // To initialize the element string with the first character (uppercase)
    std::string element;
    if (!filteredAtomName.empty()) {
        element += filteredAtomName[0];

        // If the filtered name has more than one character and the second is lowercase,
        // include the second character as part of the element (e.g., 'Ca' for Calcium)
        if (filteredAtomName.length() > 1 && std::islower(filteredAtomName[1])) {
            element += filteredAtomName[1];
        }
    }

    return element;
}

// Function to parse a single line from the PDB file

PDBAtomInfo PDBResidueParser::parsePDBLine(const std::string& line) {
    PDBAtomInfo atom;

    // To extract record type (ATOM or HETATM)
    atom.recordType = line.substr(0, 6);
    trim(atom.recordType);

    // Atom number
    atom.atomNum = std::stoi(line.substr(6, 5));

    // Atom name (e.g., "N", "CA", "HB2")
    atom.atomName = line.substr(12, 4);
    trim(atom.atomName);

    // Residue name (e.g., "MET", "SER")
    atom.resName = line.substr(17, 3);
    trim(atom.resName);

    // Chain ID (if present, otherwise use default ' ')
    char tempChainID = line[21];
    atom.chainID = (std::isalpha(tempChainID) || std::isdigit(tempChainID)) ? tempChainID : ' ';

    // Residue sequence number
    atom.resSeq = std::stoi(line.substr(22, 4));

    // Position (coordinates)
    atom.position[0] = std::stod(line.substr(30, 8)); // X
    atom.position[1] = std::stod(line.substr(38, 8)); // Y
    atom.position[2] = std::stod(line.substr(46, 8)); // Z

    // Occupancy and temperature factor
    try {
        atom.occupancy = std::stod(line.substr(54, 6));
        atom.tempFactor = std::stod(line.substr(60, 6));
    }
    catch (...) {
        // To use default values if parsing fails
        atom.occupancy = 1.0;
        atom.tempFactor = 0.0;
    }

    // Element symbol (e.g., "C", "N")
    std::string possibleElement = line.substr(76, 2);
    trim(possibleElement); // Remove trailing spaces
    possibleElement = inferElementFromAtomName(atom.atomName);

    // If element is missing or invalid, infer it from atom name
    if (possibleElement.empty() || !std::isalpha(possibleElement[0])) {
        possibleElement = inferElementFromAtomName(atom.atomName);
    }
    atom.element = possibleElement;

    return atom;
}


// Main function to parse the PDB file
void PDBResidueParser::parseFile(const string& filename, vector<Residues>& residues) {
    ifstream pdbFile(filename);
    string line;
    vector<int> atomIndices; // To store indices of atoms for the current residue
    string currentResName;
    int currentResSeq = -1;

    while (getline(pdbFile, line)) {
        if (line.substr(0, 4) == "ATOM") {
            PDBAtomInfo atom = parsePDBLine(line);
            atoms.push_back(atom);

            // To check if the residue has changed
            if (currentResSeq != atom.resSeq || currentResName != atom.resName) {
                // If there was a previous residue, process it
                if (!atomIndices.empty()) {
                    processResidue(atomIndices, currentResName, residues);
                }

                // To start a new residue
                atomIndices.clear(); // Clear the indices for the new residue
                currentResName = atom.resName;
                currentResSeq = atom.resSeq;
            }

            // Add the atom index to the current residue's indices
            atomIndices.push_back(atom.atomNum - 1);
        }
    }

    // Final residue processing
    if (!atomIndices.empty()) {
        processResidue(atomIndices, currentResName, residues);
    }

    pdbFile.close();
}


// Helper function to process each residue and identify connections
void PDBResidueParser::processResidue(vector<int> atomIndices, string currentResName, vector<Residues>& residues) {
    if (currentResName == "WAT") { // Process water molecules
        Residues wResidue;
        wResidue.resName = currentResName;
        wResidue.AllAtomsIndices = atomIndices;

        residues.push_back(wResidue);
    }
    else { // To process protein or other residues
        Residues pResidue;
        pResidue.AllAtomsIndices = atomIndices;

        pResidue.resName = currentResName;

        // To extract Hydrogen Atoms
        for (int i = pResidue.AllAtomsIndices[0]; i < pResidue.AllAtomsIndices.size(); ++i) {
            if (atoms[i].element == "H") {
                pResidue.HAtomsIndices.push_back(atoms[i].atomNum-1);
            }
            else {
                pResidue.NonHAtomsIndices.push_back(atoms[i].atomNum - 1);

            }
        }

        // Leave bond lists empty
        // Push the residue without bond information
        residues.push_back(pResidue);


    }
}
