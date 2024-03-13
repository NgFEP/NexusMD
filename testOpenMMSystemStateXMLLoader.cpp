#include "OpenMMSystemStateXMLLoader.h" // Include the state loader header
#include <iostream>
#include "OpenMM.h"

int main() {
    SystemStateXMLLoader systemstateLoader("system.xml", "state.xml"); // Load system and state

    OpenMM::System* system = systemstateLoader.getSystem();
    OpenMM::State* state = systemstateLoader.getState(); // Use a pointer to OpenMM::State

    if (system != nullptr && state != nullptr) { // Check if both system and state are successfully loaded
        std::cout << "System and State loaded successfully." << std::endl;
        // Use the system and state as needed
    }
    else {
        std::cerr << "Failed to load the system or state." << std::endl;
    }

    // Usage of system and state...

    return 0;
}
