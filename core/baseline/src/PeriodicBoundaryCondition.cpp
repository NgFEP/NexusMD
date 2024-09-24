#include "PeriodicBoundaryCondition.h"
#include <cmath>       // For abs
#include <algorithm> // For min and max
#include <stdexcept>   // For runtime_error

using namespace std;


//// Function to extract box boundaries from the initial atom positions
//void PeriodicBoundaryCondition::extractBoxBoundaries(const vector<Coords3D>& atomPositions, BoxInfo& boxInfo) {
//    if (atomPositions.empty()) {
//        throw runtime_error("No atom positions available to extract boundaries.");
//    }
//
//    boxInfo.lb = boxInfo.ub = atomPositions[0];
//    for (const auto& pos : atomPositions) {
//        for (int i = 0; i < 3; ++i) {
//            boxInfo.lb[i] = min(boxInfo.lb[i], pos[i]);
//            boxInfo.ub[i] = max(boxInfo.ub[i], pos[i]);
//        }
//    }
//
//    boxInfo.boxSize = Coords3D(abs(boxInfo.ub[0] - boxInfo.lb[0]), abs(boxInfo.ub[1] - boxInfo.lb[1]), abs(boxInfo.ub[2] - boxInfo.lb[2]));
//}


//// in case user needs to specify PBC manually
//void PeriodicBoundaryCondition::setPBC(BoxInfo& boxInfo,Coords3D& boxsize) {
//
//    boxInfo.boxSize = boxsize;
//}

void PeriodicBoundaryCondition::extractBoxBoundaries(const vector<Coords3D>& atomPositions, BoxInfo& boxInfo, const string& stateFilename) {

    bool molecule_centered = true; // to test the accuracy of the PBC I position the box in a way that molecule is at the corner of the box, otherwise if true box center is at the center of the molecule 
    bool box_crossing_molecule = false;
    bool water_availibility = true; // needs attention: do this automatically
    if (atomPositions.empty()) {
        throw runtime_error("No atom positions available to extract boundaries.");
    }
    //if (!boxInfo.boxSize[0]) {
    //    // obtain the box size
    //    boxInfo.boxSize = BaseLine::StateXMLParser::extractBoxSize(stateFilename);
    //}


    boxInfo.lb = boxInfo.ub = atomPositions[0];
    for (const auto& pos : atomPositions) {
        for (int i = 0; i < 3; ++i) {
            boxInfo.lb[i] = min(boxInfo.lb[i], pos[i]);
            boxInfo.ub[i] = max(boxInfo.ub[i], pos[i]);

        }
    }

    //updating boxInfo.boxSize to reflect the real size of the box
    boxInfo.boxSize = Coords3D(abs(boxInfo.ub[0] - boxInfo.lb[0]),
        abs(boxInfo.ub[1] - boxInfo.lb[1]),
        abs(boxInfo.ub[2] - boxInfo.lb[2]));

    if (!water_availibility) {// this is an indication of the system not having water molecules and we are dealing with a single molecule (abs(boxInfo.ub[0] - boxInfo.lb[0])) < 0.95*boxInfo.boxSize[0])
        
        if (molecule_centered){
            // Calculate the center of the molecule
            Coords3D moleculeCenter = Coords3D(
                (boxInfo.lb[0] + boxInfo.ub[0]) / 2,
                (boxInfo.lb[1] + boxInfo.ub[1]) / 2,
                (boxInfo.lb[2] + boxInfo.ub[2]) / 2
            );

            // Adjust the box so that its center matches the molecule's center
            for (int i = 0; i < 3; ++i) {
                boxInfo.lb[i] = moleculeCenter[i] - boxInfo.boxSize[i] / 2;
                boxInfo.ub[i] = moleculeCenter[i] + boxInfo.boxSize[i] / 2;
            }
        }
        else if (molecule_centered==false && box_crossing_molecule==false) {
            double marginfactor = 0.01;
            for (int i = 0; i < 3; ++i) {
                boxInfo.lb[i] -= marginfactor * boxInfo.boxSize[i];
                boxInfo.ub[i] = boxInfo.boxSize[i] + boxInfo.lb[i];// +marginfactor * boxInfo.boxSize[i]; //max(boxInfo.ub[i], pos[i]);
            }
        }
        else if (molecule_centered == false && box_crossing_molecule == true) {
            // Calculate the center of the molecule
            Coords3D moleculeCenter = Coords3D(
                (boxInfo.lb[0] + boxInfo.ub[0]) / 2,
                (boxInfo.lb[1] + boxInfo.ub[1]) / 2,
                (boxInfo.lb[2] + boxInfo.ub[2]) / 2
            );

            // Adjust the box so that its center matches the molecule's center
            for (int i = 0; i < 3; ++i) {
                boxInfo.lb[i] = moleculeCenter[i]; // lower boundary is located at the center of the molecule
                boxInfo.ub[i] = boxInfo.lb[i] + boxInfo.boxSize[i];
            }

        }
    }
    // if water molecules are in the system substracting there's no need to change the box size and .lb and .ub

    // Define a margin (e.g., 1% of the box size)
    double marginFactor = 0.01;
    Coords3D margin = Coords3D(marginFactor * boxInfo.boxSize[0],
        marginFactor * boxInfo.boxSize[1],
        marginFactor * boxInfo.boxSize[2]);

    // Adjust the boundaries by the margin
    boxInfo.lb -= margin;
    boxInfo.ub += margin;
    boxInfo.boxSize = Coords3D(abs(boxInfo.ub[0] - boxInfo.lb[0]),
        abs(boxInfo.ub[1] - boxInfo.lb[1]),
        abs(boxInfo.ub[2] - boxInfo.lb[2]));
    
} 




Coords3D PeriodicBoundaryCondition::minimumImageVector(const Coords3D& coords1, const Coords3D& coords2, const BoxInfo& boxInfo) {
    Coords3D delta = coords1 - coords2;
    for (int i = 0; i < 3; ++i) {
        if (delta[i] > 0.5 * boxInfo.boxSize[i]) {
            delta[i] -= boxInfo.boxSize[i];
        }
        else if (delta[i] < -0.5 * boxInfo.boxSize[i]) {
            delta[i] += boxInfo.boxSize[i];
        }
    }
    return delta;
}
