//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
//#include <stdexcept>   // For runtime_error
#include "pugixml.hpp"
#include "StateXMLParser.h"
// The implementation of parseSystemXML and parseStateXML goes here


using namespace std;
using namespace BaseLine;

vector<Coords3D> StateXMLParser::parseStateXML(const string& filename) {
    vector<Coords3D> atomPositions;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return atomPositions;
    }

    int id = 0;
    for (auto& position : doc.child("State").child("Positions").children("Position")) {
        Coords3D ap;
        ap[0] = position.attribute("x").as_double();
        ap[1] = position.attribute("y").as_double();
        ap[2] = position.attribute("z").as_double();
        atomPositions.push_back(ap);
    }

    return atomPositions;
}

Coords3D StateXMLParser::extractBoxSize(const string& filename) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" + filename + "] parsed with errors, Error description: " + string(result.description()) << "\n";
    }

    pugi::xml_node stateNode = doc.child("State");
    pugi::xml_node periodicBoxVectorsNode = stateNode.child("PeriodicBoxVectors");
    if (!periodicBoxVectorsNode) {
        cerr << "PeriodicBoxVectors not found in XML [" + filename + "]" << "\n";

    }

    Coords3D boxSize;
    boxSize[0] = periodicBoxVectorsNode.child("A").attribute("x").as_double();
    boxSize[1] = periodicBoxVectorsNode.child("B").attribute("y").as_double();
    boxSize[2] = periodicBoxVectorsNode.child("C").attribute("z").as_double();

    return boxSize;
}




//int main() {
//    string statefile = "state.xml";
//    auto atompositions = parsestatexml(statefile);
//
//    // example: print first atom position
//    if (!atompositions.empty()) {
//        const auto& ap = atompositions[0];
//        cout << "first atom position: id=" << ap.id << ", position=(" << ap.x << ", " << ap.y << ", " << ap.z << ")" << endl;
//    }
//
//    return 0;
//}
