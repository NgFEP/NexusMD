#include "OpenMM.h"
#include "OpenMMSystemStateXMLLoader.h"
#include <fstream>
#include <iostream>


SystemStateXMLLoader::SystemStateXMLLoader(const std::string& systemFilename, const std::string& stateFilename) {
    // Load system
    std::ifstream system_input(systemFilename.c_str());
    if (system_input) {
        systemXML = OpenMM::XmlSerializer::deserialize<OpenMM::System>(system_input);
    }
    else {
        systemXML = nullptr;
        std::cerr << "Error: Unable to open " << systemFilename << " for reading." << std::endl;
    }

    // Load state
    std::ifstream state_input(stateFilename.c_str());
    if (state_input) {
        stateXML = OpenMM::XmlSerializer::deserialize<OpenMM::State>(state_input);
    }
    else {
        stateXML = nullptr;
        std::cerr << "Error: Unable to open " << stateFilename << " for reading." << std::endl;
    }
}

SystemStateXMLLoader::~SystemStateXMLLoader() {
    delete systemXML;
    delete stateXML;
}

OpenMM::System* SystemStateXMLLoader::getSystem() const {
    return systemXML;
}

OpenMM::State* SystemStateXMLLoader::getState() const {
    return stateXML;
}
