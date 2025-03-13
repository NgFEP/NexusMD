#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

// Structure to hold constraint type and the number of constraints
struct Constraint {
    bool isHBond = false;  // True if it's an O-H bond constraint for instance
    bool isHAngle = false; // True if it's an H-O-H angle constraint for instance
    int numHBonds = 0; // Number of constraints of HBond
    int numHAngles = 0; // Number of constraints of HHAngle

};

struct ResBondInfo {
    double d;       // Ideal bond distance in nanometers
    double k;
    int bondInx;
    int p1InxLocal;
    int p2InxLocal;
    int p1InxGlobal;
    int p2InxGlobal;
    bool p1InRes = true;
    bool p2InRes = true;
};


struct CudaBondInfo {
    double d;       // Ideal bond distance in nanometers
    double k;
    int bondInx;
    int p1InxLocal;
    int p2InxLocal;
    int p1InxGlobal;
    int p2InxGlobal;
    bool waterMol = false;// Does the bond belong to water?
    bool p1InRes = true;
    bool p2InRes = true;
};
#endif // DATASTRUCTURES_H
