//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

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
