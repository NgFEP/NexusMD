// XMLParser.h
#ifndef XMLPARSER_H
#define XMLPARSER_H

#include "pugixml.hpp"
#include <vector>
#include <string>

struct TorsionParameters {
    int p1, p2, p3, p4;
    double k, phase;
    int periodicity;
};

struct AtomPosition {
    int id;
    double x, y, z;
};

class XMLParser {
public:
    // Parses the system.xml file and returns a vector of TorsionForce structures
    static std::vector<TorsionParameters> parseSystemXML(const std::string& filename);

    // Parses the state.xml file and returns a vector of AtomPosition structures
    static std::vector<AtomPosition> parseStateXML(const std::string& filename);
};

#endif // XMLPARSER_H
