// SystemXMLParser.h
#ifndef SYSTEMXMLPARSER_H
#define SYSTEMXMLPARSER_H

#include "pugixml.hpp"
#include <vector>
#include <string>
using namespace std;
namespace BaseLine {

    struct TorsionParameters {
        int p1, p2, p3, p4;
        double k, phase;
        int periodicity;
    };

    class SystemXMLParser {
    public:

        static vector<double> MassesParser(const string& filename);

        // Parses the state.xml file and returns a vector of AtomPosition structures
        static vector<TorsionParameters> PeriodicTorsionParser(const string& filename);
    };
}

#endif // SYSTEMXMLPARSER_H
