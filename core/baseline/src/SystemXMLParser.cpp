//#include "stdafx.h"

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "pugixml.hpp"
#include "SystemXMLParser.h"
// The implementation of parseSystemXML and parseStateXML goes here


using namespace std;
using namespace BaseLine;


vector<double> SystemXMLParser::MassesParser(const string& filename) {
    vector<double> Masses;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return Masses;
    }
    //Navigating <Particle> elements under <Particles>
    for (auto& particle : doc.child("System").child("Particles").children("Particle")) {
        double mass = particle.attribute("mass").as_double();
        Masses.push_back(mass);
    }

    return Masses;

}

vector<PTorsionParams> SystemXMLParser::PTorsionParser(const string& filename) {
    vector<PTorsionParams> torsionParameters;
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
                PTorsionParams tp;
                tp.p1 = torsion.attribute("p1").as_int();
                tp.p2 = torsion.attribute("p2").as_int();
                tp.p3 = torsion.attribute("p3").as_int();
                tp.p4 = torsion.attribute("p4").as_int();
                tp.k = torsion.attribute("k").as_double();// k type is double
                tp.periodicity = torsion.attribute("periodicity").as_int();// periodicity type is int
                tp.phase = torsion.attribute("phase").as_double();// phase type is double
                torsionParameters.push_back(tp);
            }
        }
    }

    return torsionParameters;
}

vector<HBondParams> SystemXMLParser::HBondParser(const string& filename) {
    vector<HBondParams> bondParameters;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return bondParameters;
    }

    // Navigate to the <Force> elements and check for the type attribute to be "HarmonicBondForce"
    for (auto& force : doc.child("System").child("Forces").children("Force")) {
        string forceType = force.attribute("type").as_string();
        if (forceType == "HarmonicBondForce") {
            // Once the correct <Force> element is found, iterate through its <Bond> children
            for (auto& bond : force.child("Bonds").children("Bond")) {
                HBondParams bp;
                bp.p1 = bond.attribute("p1").as_int();
                bp.p2 = bond.attribute("p2").as_int();
                bp.d = bond.attribute("d").as_double();
                bp.k = bond.attribute("k").as_double();
                bondParameters.push_back(bp);
            }
        }
    }

    return bondParameters;
}


vector<HAngleParams> SystemXMLParser::HAngleParser(const string& filename) {
    vector<HAngleParams> angleParameters;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());

    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return angleParameters;
    }

    // Navigate to the <Force> elements and check for the type attribute to be "HarmonicAngleForce"
    for (auto& force : doc.child("System").child("Forces").children("Force")) {
        string forceType = force.attribute("type").as_string();
        if (forceType == "HarmonicAngleForce") {
            // Once the correct <Force> element is found, iterate through its <Angle> children
            for (auto& angle : force.child("Angles").children("Angle")) {
                HAngleParams ap;
                ap.p1 = angle.attribute("p1").as_int();
                ap.p2 = angle.attribute("p2").as_int();
                ap.p3 = angle.attribute("p3").as_int();
                ap.a = angle.attribute("a").as_double();
                ap.k = angle.attribute("k").as_double();
                angleParameters.push_back(ap);
            }
        }
    }

    return angleParameters;
}


NonbondedParams SystemXMLParser::NonBondedParser(const string& filename) {
    NonbondedParams params;
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if (!result) {
        cerr << "XML [" << filename << "] parsed with errors, Error description: " << result.description() << "\n";
        return params;
    }

    auto forceNode = doc.child("System").find_child_by_attribute("Force", "type", "NonbondedForce");
    if (!forceNode) {
        cerr << "No NonbondedForce found in the XML.\n";
        return params;
    }
    else {
        // Existing code for parsing attributes
        params.alpha = forceNode.attribute("alpha").as_double();
        params.cutoff = forceNode.attribute("cutoff").as_double();
        params.dispersionCorrection = forceNode.attribute("dispersionCorrection").as_int();
        params.ewaldTolerance = forceNode.attribute("ewaldTolerance").as_double();
        params.exceptionsUsePeriodic = forceNode.attribute("exceptionsUsePeriodic").as_int();
        params.forceGroup = forceNode.attribute("forceGroup").as_int();
        params.includeDirectSpace = forceNode.attribute("includeDirectSpace").as_int();
        params.ljAlpha = forceNode.attribute("ljAlpha").as_int();
        params.ljnx = forceNode.attribute("ljnx").as_int();
        params.ljny = forceNode.attribute("ljny").as_int();
        params.ljnz = forceNode.attribute("ljnz").as_int();
        params.method = forceNode.attribute("method").as_int();
        params.nx = forceNode.attribute("nx").as_int();
        params.ny = forceNode.attribute("ny").as_int();
        params.nz = forceNode.attribute("nz").as_int();
        params.recipForceGroup = forceNode.attribute("recipForceGroup").as_int();
        params.rfDielectric = forceNode.attribute("rfDielectric").as_double();
        params.switchingDistance = forceNode.attribute("switchingDistance").as_double();
        params.useSwitchingFunction = forceNode.attribute("useSwitchingFunction").as_int();
        params.version = forceNode.attribute("version").as_int();

        // Parse particles
        for (auto& particle : forceNode.child("Particles").children("Particle")) {
            NonbondedParticle p;
            p.q = particle.attribute("q").as_double();
            p.sig = particle.attribute("sig").as_double();
            p.eps = particle.attribute("eps").as_double();
            params.particles.push_back(p);
        }

        // Parse exceptions
        for (auto& exception : forceNode.child("Exceptions").children("Exception")) {
            NonbondedException e;
            e.p1 = exception.attribute("p1").as_int();
            e.p2 = exception.attribute("p2").as_int();
            e.q = exception.attribute("q").as_double();
            e.sig = exception.attribute("sig").as_double();
            e.eps = exception.attribute("eps").as_double();
            params.exceptions.push_back(e);
        }

        return params;
    }


}




//
//int main() {
//    string systemFile = "system.xml";
//    auto torsionForces = PeriodicTorsionParser(systemFile);
//    auto masses = MassesParser(systemFile);
//
//    // Example: print first torsion force parameters
//    if (!torsionForces.empty()) {
//        const auto& tf = torsionForces[0];
//        cout << "First Torsion Force: p1=" << tf.p1 << ", p2=" << tf.p2
//            << ", p3=" << tf.p3 << ", p4=" << tf.p4 << ", k=" << tf.k
//            << ", periodicity=" << tf.periodicity << ", phase=" << tf.phase << endl;
//    }
//    else {
//        cout << "system.xml is empty!!!" << endl;
//    }
//
//    //printing the first particle's mass
//    if (!masses.empty()) {
//        cout << "First Particle Mass:" << masses[0] << endl;
//    }
//    else {
//        cout << "No Particle mass was found in " << systemFile << endl;
//    }
//
//
//    return 0;
//}
