// here the resources and processors are allocated based on the User selections to choose among a CPU or CUDA implementation
#include "TaskDispatcher.h"
#include "BaseLineTaskDispatcher.h"
#include "CUDATaskDispatcher.h"

std::unique_ptr<TaskDispatcher> TaskDispatcher::CreateDispatcher(const std::string& processorChoice) {
    if (processorChoice == "CPU") {
        return std::make_unique<BaseLineTaskDispatcher>();
    }
    else if (processorChoice == "CUDA") {
        return std::make_unique<CudaTaskDispatcher>();
    }
    else {
        std::cerr << "Invalid processor type. Defaulting to CPU." << std::endl;
        return std::make_unique<BaseLineTaskDispatcher>();
    }
}
