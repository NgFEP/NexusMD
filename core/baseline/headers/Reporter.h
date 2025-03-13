#ifndef REPORTER_H
#define REPORTER_H

#include "Coords3D.h"
#include <vector>
#include <string>
#include "SystemXMLParser.h"
#include <set>

namespace BaseLine {

    struct PDBAtom {
        int atomNum;
        std::string atomName;
        std::string resName;
        char chainID;
        int resSeq;
        Coords3D coords;
        double occupancy;
        double tempFactor;
        std::string element;
        std::string recordType; // distinguishes ATOM from HETATM
    };

    class Reporter {
    public:
        Reporter();
        // Update report function to no longer require clearFile parameter
        // for test report and shows positions, velocities and forces of all atoms and for all the steps, used for test purposes
        void TestPVFReport(const std::string& filename, const std::vector<Coords3D>& positions,
            const std::vector<Coords3D>& velocities, const std::vector<Coords3D>& forces, int step, const std::vector<PTorsionParams>& torsionParams, const std::vector<BondParams>& bondParams, const std::vector<AngleParams>& angleParams);
        void TotalEnergyReport(const std::string& baseFilename, double& totalKEnergy, double& totalPEnergy, double& totalEnergy, int step);

        // Main function to read, process, and write PDB file
        void pdbOutputGeneratorPart1(const std::string& inputFilename, const std::string& outputFilename, std::vector<PDBAtom>& outputTemplate);
        void pdbOutputGeneratorPart2(const std::string& outputFilename, std::vector<PDBAtom>& outputTemplate, const std::vector<Coords3D>& positions, int ModelNum);
        void EnsembleReport(const std::string& baseFilename, double& totalKEnergy, double& calcTemp, double& volume, double& density, int& step);


    private:
        std::string _PBCLine = "";
        // To set containing standard amino acid residues
        const std::set<std::string> standardAminoAcids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        };
        const std::set<std::string> dnaNucleotides = {
            "DA", // Adenine
            "DT", // Thymine
            "DG", // Guanine
            "DC"  // Cytosine
        };
        const std::set<std::string> rnaNucleotides = {
            "A", // Adenine
            "U", // Uracil
            "G", // Guanine
            "C"  // Cytosine
        };

        const std::set<std::string> standardResidues = {
            // Standard Amino Acids
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",

            // DNA Nucleotides
            "DA", "DT", "DG", "DC",

            // RNA Nucleotides
            "A", "U", "G", "C"
        };

        //static ofstream outputFile;
        std::string getCurrentDate();

        // Extracts the element from the atom name
        std::string extractElement(const std::string& atomName);

        PDBAtom parsePDBLine(const std::string& line);

        void outputTemplateExtractor(const std::string& inputFilename, std::vector<PDBAtom>& outputTemplate);

        // Helper function to format a PDB line
        std::string formatPDBLine(const PDBAtom& atom);

    };
}

#endif // REPORTER_H
