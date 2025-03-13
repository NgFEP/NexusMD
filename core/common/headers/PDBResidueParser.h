#ifndef PDBRESIDUEPARSER_H
#define PDBRESIDUEPARSER_H

#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Coords3D.h"
#include <iomanip>
#include <cctype>
#include <algorithm>
#include <optional>
#include "DataStructures.h"
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


struct Molecules {
    std::vector<int> AllAtomsIndices; // List of atom indices
    std::string type; // Polymer (protein) and Non-polymer (ligand, water) 
};


// Unified host Residue Struct for storing residue information of polymer and non-polymer 
struct Residues {
    std::vector<int> AllAtomsIndices; // List of atom indices
    std::string resName;
    std::vector<int> HAtomsIndices; // List of Hydrogen atom indices
    std::vector<int> NonHAtomsIndices; // List of Non Hydrogen atom indices
    std::optional<std::vector<ResBondInfo>> AllBondsIndices; // List of all bond indices
    //std::optional<std::vector<int>> AllBondsIndices; // List of all bond indices
    std::optional<std::vector<ResBondInfo>> HBondsIndices; // List of Hydrogen bond indices
    std::optional<std::vector<ResBondInfo>> NonHBondsIndices;// List of non-Hydrogen bond indices
    std::optional<std::vector<int>> NonResBondAtoms; // List of all atoms in bonds that not belong to the current residue
    std::optional<std::vector<int>> AngleIndices;
    std::optional<std::vector<int>> TorsionIndices;
    int resMemSize = 0; // Residue memory size
};



//// Struct for storing connection between residues
//struct Connection {
//    int atomC; // C atom number of previous residue
//    int atomN; // N atom number of current residue
//    //int Bond = std::nullopt;  // Bond index between residues
//    std::optional<int> BondIndex; // bond index, each connection has a single bond
//};

//// Struct for storing water information
//struct WResidues {
//    int lowBound;
//    int highBound;
//    //std::vector<int> AtomsIDs; // List of atom indices
//    //std::vector<int> HAtomsIDs; // List of Hydrogen atom indices in water
//    //std::vector<int> NonHAtomsIDs; // List of Non Hydrogen atom indices
//    //std::unordered_set<int> HAtomsIDs;
//    std::optional<std::vector<int>> BondsIndices; // List of Hydrogen bond indices
//};

// Main PDB Parser class
class PDBResidueParser {
public:
    PDBResidueParser();
    ~PDBResidueParser();
    void parseFile(const std::string& filename, std::vector<Residues>& Residues, std::vector<Molecules>& molecules);
    //const std::vector<PResidues>& getResidues() const { return residues; }
    //const std::vector<Connection>& getConnections() const { return connections; }
    //const std::vector<wResidues>& getWaters() const { return wResidues; }

private:
    //std::string filename;
    std::vector<PDBAtomInfo> atoms;
    //std::vector<PResidues> residues;
    //std::vector<Connection> connections;
    //std::vector<wResidues> wResidues;

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
