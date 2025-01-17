#ifndef PERIODICBOUNDARYCONDITION_H
#define PERIODICBOUNDARYCONDITION_H

#include "Coords3D.h"
#include "vector"
#include "StateXMLParser.h"


class PeriodicBoundaryCondition {
public:
    struct BoxInfo {
        Coords3D lb; // lower boundary
        Coords3D ub; // upper boundary
        Coords3D boxSize;

        BoxInfo() : lb(0.0, 0.0, 0.0), ub(0.0, 0.0, 0.0), boxSize(0.0, 0.0, 0.0) {}
        BoxInfo(const Coords3D& lower, const Coords3D& upper)
            : lb(lower), ub(upper), boxSize(std::abs(upper[0] - lower[0]), std::abs(upper[1] - lower[1]), std::abs(upper[2] - lower[2])) {}
    };
        
    static void extractBoxBoundaries(const std::vector<Coords3D>& atomPositions, BoxInfo& boxInfo, const std::string& stateFilename);
    static Coords3D minimumImageVector(const Coords3D& coords1, const Coords3D& coords2, const BoxInfo& boxInfo);

};


#endif // PERIODICBOUNDARYCONDITION_H
