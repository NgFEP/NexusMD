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

    struct HAngleParams {// HAngleParams HarmonicAngleParams
        int p1, p2, p3; // Indices of the particles forming the angle
        double a; // Ideal angle in radians
        double k; // Force constant
    };


    class SystemXMLParser {
    public:

        static std::vector<double> MassesParser(const std::string& filename);

        // Parses the state.xml file and returns a vector of PTorsionParams structure
        static std::vector<PTorsionParams> PTorsionParser(const std::string& filename);

        // Parses the state.xml file and returns a vector of HAngleParams structure
        static std::vector<HAngleParams> HAngleParser(const std::string& filename);

    };

}

#endif // SYSTEMXMLPARSER_H
