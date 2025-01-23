#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

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
    bool waterMol = false;// does the bond belong to water?
    bool p1InRes = true;
    bool p2InRes = true;
};
#endif // DATASTRUCTURES_H
