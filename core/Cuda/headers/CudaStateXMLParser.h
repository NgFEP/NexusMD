// StateXMLParser.h
#ifndef CUDASTATEXMLPARSER_H
#define CUDASTATEXMLPARSER_H

#include "pugixml.hpp"
#include <vector>
#include <string>
#include "Coords3D.h"

namespace Cuda {


    class StateXMLParser {
    public:

        // Parses the state.xml file and returns a vector of AtomPosition structures
        static std::vector<Coords3D> parseStateXML(const std::string& filename);
        static Coords3D extractBoxSize(const std::string& filename);

    };
}

#endif // CUDASTATEXMLPARSER_H
