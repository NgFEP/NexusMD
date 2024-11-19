#ifndef PDBSEGMENTPARSER_H
#define PDBSEGMENTPARSER_H

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

// Struct for storing residue information
struct PResidues {
    int lowBound;
    int highBound;
    std::string resName;
    std::vector<int> HAtomsIDs; // List of Hydrogen atom indices
    //std::unordered_set<int> HAtomsIDs;
    std::optional<std::vector<int>> AllBondsIndices; // List of all bond indices
    std::optional<std::vector<int>> HBondsIndices; // List of Hydrogen bond indices
    std::optional<std::vector<int>> NonHBondsIndices;// List of non-Hydrogen bond indices
    std::optional<std::vector<int>> AngleIndices;
    std::optional<std::vector<int>> TorsionIndices;
};

// Struct for storing connection between residues
struct Connection {
    int atomC; // C atom number of previous residue
    int atomN; // N atom number of current residue
    //int Bond = std::nullopt;  // Bond index between residues
    std::optional<int> BondIndex; // bond index, each connection has a single bond
};

// Struct for storing water information
struct WaterMols {
    int lowBound;
    int highBound;
    std::vector<int> HAtomsIDs; // List of Hydrogen atom indices in water
    //std::unordered_set<int> HAtomsIDs;
    std::optional<std::vector<int>> BondsIndices; // List of Hydrogen bond indices
};

// Main PDB Parser class
class PDBSegmentParser {
public:
    PDBSegmentParser();
    ~PDBSegmentParser();
    void parseFile(const std::string& filename, std::vector<PResidues>& residues, std::vector<Connection>& connections, std::vector<WaterMols>& waterMols);
    //const std::vector<PResidues>& getResidues() const { return residues; }
    //const std::vector<Connection>& getConnections() const { return connections; }
    //const std::vector<WaterMols>& getWaters() const { return waterMols; }

private:
    //std::string filename;
    std::vector<PDBAtomInfo> atoms;
    //std::vector<PResidues> residues;
    //std::vector<Connection> connections;
    //std::vector<WaterMols> waterMols;

    const std::set<std::string> standardAminoAcids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    PDBAtomInfo parsePDBLine(const std::string& line);
    void processResidue(int startIdx, int endIdx, const std::string& currentResName, std::vector<PResidues>& residues, std::vector<Connection>& connections, std::vector<WaterMols>& waterMols);

    // Utility function to infer element from atom name
    std::string inferElementFromAtomName(const std::string& atomName);
    // Utility function to trim whitespace
    void trim(std::string& str);

};

#endif // PDBSEGMENTPARSER_H
