#ifndef OPENMM_SYSTEM_STATE_XML_LOADER_H
#define OPENMM_SYSTEM_STATE_XML_LOADER_H

#include "OpenMM.h"
#include <string>

class SystemStateXMLLoader {
public:
    SystemStateXMLLoader(const std::string& systemFilename, const std::string& stateFilename); // Constructor
    ~SystemStateXMLLoader(); // Destructor

    OpenMM::System* getSystem() const; // Accessor for the system pointer
    OpenMM::State* getState() const; // Accessor for the state pointer

private:
    OpenMM::System* systemXML; // Pointer to the System object
    OpenMM::State* stateXML; // Pointer to the State object
};

#endif // OPENMM_SYSTEM_STATE_XML_LOADER_H
