#include "OpenMMSystemXMLLoader.h"
#include <iostream>

int main() {
    SystemXMLLoader loader("system.xml");// loader as an object for SystemXMLLoader class
    OpenMM::System* system = loader.getSystem();

    if (system != nullptr) {
        // Use the system as needed
        std::cout << "System loaded successfully." << std::endl;
    }
    else {
        std::cerr << "Failed to load the system." << std::endl;
    }

    //here is the rest of the usage of systemXML from OpenMM for the NexaBind simulator engine







    //

    // No need to manually delete the system; the SystemXMLLoader destructor handles it
    return 0;
}
