#include "TaskDispatcher.h"
#include "BaseLineTaskDispatcher.h"
#include "CUDATaskDispatcher.h"


// Factory method to create the appropriate dispatcher
std::shared_ptr<TaskDispatcher> TaskDispatcher::CreateDispatcher(const std::string& processorChoice) {
    if (processorChoice == "CPU") {
        return std::make_shared<BaseLineTaskDispatcher>();
    }
    else if (processorChoice == "CUDA") {
        return std::make_shared<CudaTaskDispatcher>();
    }
    else {
        std::cerr << "Invalid processor type. Defaulting to CPU." << std::endl;
        return std::make_shared<BaseLineTaskDispatcher>();
    }
}

// Initialize default simulation parameters
void TaskDispatcher::initializeParameters() {
    // Initialize atom positions
    atomPositions = {};

    // Initialize particle masses
    masses = {};

    // Initialize torsion, bond, and angle parameters
    torsionParams = {};
    bondParams = {};
    angleParams = {};

    // Initialize nonbonded parameters
    nonbondedParams = {};

    // Initialize periodic boundary condition box
    boxInfo.boxSize = {};
}
