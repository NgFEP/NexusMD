//#ifndef REPORTER_H
//#define REPORTER_H
//
//#include "Coords3D.h"
//#include <vector>
//#include <string>
//
//namespace Cuda {
//
//    class Reporter {
//    public:
//        // Constructs the Reporter object.
//        Reporter() = default;
//
//        // Reports simulation data.
//        // The clearFile parameter determines whether to clear existing file data (true for the first step).
//        static void report(const string& filename, const vector<Coords3D>& positions,
//            const vector<Coords3D>& velocities, const vector<Coords3D>& forces,
//            int step, bool clearFile);
//    };
//
//} // namespace Cuda
//
//#endif // REPORTER_H


#ifndef CUDAREPORTER_H
#define CUDAREPORTER_H

#include "Coords3D.h"
#include <vector>
#include <string>
#include "SystemXMLParser.h"
#include <set>

namespace Cuda {

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
        std::string recordType; // Added to distinguish ATOM from HETATM
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
        //void pdbOutputGenerator(const std::string& inputFilename, const std::string& outputFilename, std::vector<PDBAtom>& outputTemplate,const std::vector<Coords3D>& positions, int step);
        void pdbOutputGeneratorPart1(const std::string& inputFilename, const std::string& outputFilename, std::vector<PDBAtom>& outputTemplate);
        void pdbOutputGeneratorPart2(const std::string& outputFilename, std::vector<PDBAtom>& outputTemplate, const std::vector<Coords3D>& positions, int ModelNum);


    private:
        std::string _PBCLine = "";
        // Set containing standard amino acid residues
        const std::set<std::string> standardAminoAcids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
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

#endif // CUDAREPORTER_H
