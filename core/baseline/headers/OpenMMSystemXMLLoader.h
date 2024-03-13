#ifndef OpenMM_SYSTEM_XML_LOADER_H
#define OpenMM_SYSTEM_XML_LOADER_H

#include "OpenMM.h"
#include <string>

class SystemXMLLoader {
public:
    SystemXMLLoader(const std::string& filename); // Constructor
    ~SystemXMLLoader(); // Destructor

    OpenMM::System* getSystem() const; // Accessor for the system pointer
    //marking getSystem() as a const member function is a best practice for accessor functions, ensuring that the function can be used in a wide variety of contexts without risking unintended modifications to the object's state.
private:
    OpenMM::System* systemXML; // Pointer to the System object
};

#endif // OpenMM_SYSTEM_XML_LOADER_H
