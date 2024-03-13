#include "OpenMM.h"
#include "OpenMMSystemXMLLoader.h"
#include <fstream>
#include <iostream>

SystemXMLLoader::SystemXMLLoader(const std::string& filename) {// class SystemXMLLoader and constructor :: SystemXMLLoader
    std::ifstream system_input(filename.c_str());
    if (system_input) {
        systemXML = OpenMM::XmlSerializer::deserialize<OpenMM::System>(system_input);
    }
    else {
        systemXML = nullptr;
        // Handle the error as appropriate
        std::cerr << "Error: Unable to open " << filename << " for reading." << std::endl;
    }
}

SystemXMLLoader::~SystemXMLLoader() {
    delete systemXML;
}

OpenMM::System* SystemXMLLoader::getSystem() const {
    return systemXML;
}
