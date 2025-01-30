#ifndef PDBRESIDUEPARSER_H
#define PDBRESIDUEPARSER_H

#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Coords3D.h"
#include <optional>
#include <iomanip>
#include <cctype>
#include <algorithm>
//#include <unordered_set>

// Struct for storing atom information
struct PDBAtomInfo {
    std::string recordType; // e.g., ATOM or HETATM
    int atomNum;
    std::string atomName;
    std::string resName;
    char chainID;
    int resSeq;
    Coords3D position; // Using Coords3D for coordinates
    double occupancy;
    double tempFactor;
    std::string element;
};

// Unified host Residue Struct for storing residue information of Protein, water, ...
struct Residues {
    std::vector<int> AllAtomsIndices; // List of atom indices
    std::string resName;
    std::vector<int> HAtomsIndices; // List of Hydrogen atom indices
    std::vector<int> NonHAtomsIndices; // List of Non Hydrogen atom indices
    std::optional<std::vector<int>> AllBondsIndices; // List of all bond indices
    std::optional<std::vector<int>> HBondsIndices; // List of Hydrogen bond indices
    std::optional<std::vector<int>> NonHBondsIndices;// List of non-Hydrogen bond indices
    std::optional<std::vector<int>> NonResBondAtoms; // List of all atoms in bonds that not belong to the current residue
    std::optional<std::vector<int>> AngleIndices;
    std::optional<std::vector<int>> TorsionIndices;
    int resMemSize = 0; // residue memory size
};



// Main PDB Parser class
class PDBResidueParser {
public:
    PDBResidueParser();
    ~PDBResidueParser();
    void parseFile(const std::string& filename, std::vector<Residues>& Residues);

private:
    std::vector<PDBAtomInfo> atoms;

    const std::set<std::string> standardAminoAcids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    PDBAtomInfo parsePDBLine(const std::string& line);
    void processResidue(std::vector<int> atomIndices, std::string currentResName, std::vector<Residues>& Residues);

    // Utility function to infer element from atom name
    std::string inferElementFromAtomName(const std::string& atomName);
    // Utility function to trim whitespace
    void trim(std::string& str);

};

#endif // PDBRESIDUEPARSER_H
