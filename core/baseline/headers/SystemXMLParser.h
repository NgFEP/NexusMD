// SystemXMLParser.h
#ifndef SYSTEMXMLPARSER_H
#define SYSTEMXMLPARSER_H

#include "pugixml.hpp"
#include <vector>
#include <string>

namespace BaseLine {

    struct PTorsionParams {// used for Periodic Torsion Force (PTF)
        int p1, p2, p3, p4;
        double k, phase;
        int periodicity;
    };

    struct HBondParams {
        int p1, p2;     // Indices of the two particles forming the bond
        double d;       // Ideal bond distance in nanometers
        double k;       // Force constant in kJ/mol/nm^2
    };

    struct HAngleParams {// HAngleParams HarmonicAngleParams
        int p1, p2, p3; // Indices of the particles forming the angle
        double a; // Ideal angle in radians
        double k; // Force constant
    };


    struct NonbondedParticle {
        double q, sig, eps;
    };

    struct NonbondedException {
        int p1, p2;
        double q, sig, eps;
    };

    struct NonbondedParams {
        double alpha, cutoff, ewaldTolerance, rfDielectric, switchingDistance;
        int dispersionCorrection, exceptionsUsePeriodic, forceGroup, includeDirectSpace, ljAlpha, ljnx, ljny, ljnz, method, nx, ny, nz, recipForceGroup, useSwitchingFunction, version;
        std::vector<NonbondedParticle> particles;
        std::vector<NonbondedException> exceptions;
    };




    class SystemXMLParser {
    public:

        static std::vector<double> MassesParser(const std::string& filename);

        // Parses the system.xml file and returns a vector of PTorsionParams structure
        static std::vector<PTorsionParams> PTorsionParser(const std::string& filename);

        // Parses the system.xml file and returns a vector of HBondParams structure
        static std::vector<HBondParams> HBondParser(const std::string& filename);

        // Parses the system.xml file and returns a vector of HAngleParams structure
        static std::vector<HAngleParams> HAngleParser(const std::string& filename);

        static NonbondedParams NonBondedParser(const std::string& filename);  // New function



    };

}

#endif // SYSTEMXMLPARSER_H
