#include "Reporter.h"
#include <iostream>
#include <fstream>
#include <iomanip> // For std::fixed, std::setprecision, and std::setw

using namespace BaseLine;
using namespace std;

void Reporter::report(const std::string& filename, const std::vector<Coords3D>& positions, const std::vector<Coords3D>& velocities, const std::vector<Coords3D>& forces, int step, bool clearFile) {
    // Determine file open mode based on whether clearing the file or appending
    ios_base::openmode fileMode = clearFile ? (ios::out) : (ios::out | ios::app);

    // Open the file with the selected mode
    ofstream outputFile(filename, fileMode);

    if (!outputFile.is_open()) {
        cerr << "Unable to open the file: " << filename << endl;
        return;
    }

    // Write headers only if clearing the file
    if (clearFile) {
        outputFile << "Step\tID\tPosition\t\t\t\t\tVelocity\t\t\t\t\tForce\n";
    }

    // Set fixed number of decimal places for floating-point output
    outputFile << fixed << setprecision(6);

    int numAtoms = min(10, static_cast<int>(positions.size()));
    for (int i = 0; i < numAtoms; ++i) {
        outputFile << step << "\t" << i + 1
            << "\t(" << setw(10) << positions[i][0] << ", " << setw(10) << positions[i][1] << ", " << setw(10) << positions[i][2] << ")\t"
            << "(" << setw(10) << velocities[i][0] << ", " << setw(10) << velocities[i][1] << ", " << setw(10) << velocities[i][2] << ")\t"
            << "(" << setw(10) << forces[i][0] << ", " << setw(10) << forces[i][1] << ", " << setw(10) << forces[i][2] << ")\n";
    }

    outputFile.close();
}
