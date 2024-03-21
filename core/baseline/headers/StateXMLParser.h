// StateXMLParser.h
#ifndef STATEXMLPARSER_H
#define STATEXMLPARSER_H

#include "pugixml.hpp"
#include <vector>
#include <string>
#include "Coords3D.h"

namespace BaseLine {


    class StateXMLParser {
    public:

        // Parses the state.xml file and returns a vector of AtomPosition structures
        static std::vector<Coords3D> parseStateXML(const std::string& filename);
    };
}

#endif // STATEXMLPARSER_H
