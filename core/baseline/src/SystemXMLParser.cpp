#include "stdafx.h"
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
