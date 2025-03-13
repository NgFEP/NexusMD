#pragma once

#include "TaskDispatcher.h"
#include "BaseLineTaskDispatcher.h"
#include "CudaTaskDispatcher.h"

std::shared_ptr<TaskDispatcher> TaskDispatcher::CreateDispatcher(const std::string& processorChoice) {
    if (processorChoice == "CPU" || processorChoice == "Cpu" || processorChoice == "cpu") {
        return std::make_shared<BaseLineTaskDispatcher>();
    }
    else if (processorChoice == "CUDA" || processorChoice == "Cuda" || processorChoice == "cuda") {
        return std::make_shared<CudaTaskDispatcher>();
    }
    else {
        std::cerr << "Invalid processor type. Defaulting to CPU." << std::endl;
        return std::make_shared<BaseLineTaskDispatcher>();
    }
}

void TaskDispatcher::initializeParameters() {
    atomPositions = {};
    masses = {};
    torsionParams = {};
    bondParams = {};
    angleParams = {};
    nonbondedParams = {};
}
