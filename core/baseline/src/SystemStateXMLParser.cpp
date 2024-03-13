#include "pugixml.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "SystemStateXMLParser.h"
// The implementation of parseSystemXML and parseStateXML goes here


using namespace std;


vector<TorsionParameters> parseSystemXML(const string& filename) {
    vector<TorsionParameters> torsionParameters;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return torsionParameters;
    }

    // Navigate to the <Force> elements and check for the type attribute to be "PeriodicTorsionForce"
    for (auto& force : doc.child("System").child("Forces").children("Force")) {
        string forceType = force.attribute("type").as_string();
        if (forceType == "PeriodicTorsionForce") {
            // Once the correct <Force> element is found, iterate through its <Torsion> children
            for (auto& torsion : force.child("Torsions").children("Torsion")) {
                TorsionParameters tp;
                tp.p1 = torsion.attribute("p1").as_int();
                tp.p2 = torsion.attribute("p2").as_int();
                tp.p3 = torsion.attribute("p3").as_int();
                tp.p4 = torsion.attribute("p4").as_int();
                tp.k = torsion.attribute("k").as_double();
                tp.periodicity = torsion.attribute("periodicity").as_int();
                tp.phase = torsion.attribute("phase").as_double();
                torsionParameters.push_back(tp);
            }
        }
    }

    return torsionParameters;
}


vector<AtomPosition> parseStateXML(const string& filename) {
    vector<AtomPosition> atomPositions;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return atomPositions;
    }

    int id = 0;
    for (auto& position : doc.child("State").child("Positions").children("Position")) {
        AtomPosition ap;
        ap.id = id++;
        ap.x = position.attribute("x").as_double();
        ap.y = position.attribute("y").as_double();
        ap.z = position.attribute("z").as_double();
        atomPositions.push_back(ap);
    }

    return atomPositions;
}

int main() {
    string systemFile = "system.xml";
    string stateFile = "state.xml";
    auto torsionForces = parseSystemXML(systemFile);
    auto atomPositions = parseStateXML(stateFile);

    // Example: print first torsion force parameters
    if (!torsionForces.empty()) {
        const auto& tf = torsionForces[0];
        cout << "First Torsion Force: p1=" << tf.p1 << ", p2=" << tf.p2
            << ", p3=" << tf.p3 << ", p4=" << tf.p4 << ", k=" << tf.k
            << ", periodicity=" << tf.periodicity << ", phase=" << tf.phase << endl;
    }
    else {
        cout << "system.xml is empty!!!" << endl;
    }

    // Example: print first atom position
    if (!atomPositions.empty()) {
        const auto& ap = atomPositions[0];
        cout << "First Atom Position: id=" << ap.id << ", position=(" << ap.x << ", " << ap.y << ", " << ap.z << ")" << endl;
    }

    return 0;
}
