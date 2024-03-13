#include "PeriodicTorsionForce.h"
#include "Coords3D.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "SystemStateXMLParser.h" //because of its Torsionparameters

#ifndef PI
#define PI 3.14159265358979323846
#endif

using namespace std;
using namespace BaseLine; // Coords3D is within the BaseLine namespace

PeriodicTorsionForce::PeriodicTorsionForce() {
}

PeriodicTorsionForce::~PeriodicTorsionForce() {
}

// Calculate distance between two 3-d cartesian coordinates
double PeriodicTorsionForce::atom_distance(const Coords3D& coords1, const Coords3D& coords2) {
    Coords3D diff = coords2 - coords1; // Utilizing Coords3D operator overloading
    return sqrt(diff.dot(diff)); // Using the dot method from Coords3D
}


// Calculate unit cross product between two vectors in each plane to find the perpendicular normal vector to the plane
Coords3D PeriodicTorsionForce::plane_norm_vec(const Coords3D& uvec1, const Coords3D& uvec2) {
    double cos_12 = uvec1.dot(uvec2);
    cos_12 = max(min(cos_12, 1.0), -1.0);//making sure cos is between -1 and +1
    double sin_12 = sqrt(1 - cos_12 * cos_12);
    Coords3D ucp = uvec1.cross(uvec2) / sin_12;

    return ucp;
}

// Calculate torsion angle between four 3-d cartesian coordinates, we need torsion angle in rad to later on 
//calculate periodic torsion force and energy
double PeriodicTorsionForce::torsion_angle(const Coords3D& coords1, const Coords3D& coords2,
    const Coords3D& coords3, const Coords3D& coords4) {
    Coords3D u21 = (coords2 - coords1)/ atom_distance(coords2, coords1);
    Coords3D u23 = (coords2 - coords3)/ atom_distance(coords2, coords3);
    Coords3D u32 = (coords3 - coords2)/ atom_distance(coords3, coords2);
    Coords3D u34 = (coords3 - coords4)/ atom_distance(coords3, coords4);

    Coords3D u21c23 = plane_norm_vec(u21, u23);
    Coords3D u32c34 = plane_norm_vec(u32, u34);

    double dp = u21c23.dot(u32c34);//dp=dot product of two planes normal vectors
    dp = max(min(dp, 1.0), -1.0);//making sure dp which is the cos of torsion angle is between -1 and +1
    double sign = (u21c23.dot(u34) >= 0) ? 1.0 : -1.0;
    double t1234 = sign * acos(dp); // t1234 is now in radians


    // Calculate vectors between atoms
    Coords3D v12 = u21;
    Coords3D v23 = u32;
    Coords3D v34 = -u34;

    // Calculate normals to the planes defined by first three and last three atoms
    Coords3D n1 = v12.cross(v23).normalize();
    Coords3D n2 = v23.cross(v34).normalize();

    // Calculate vector perpendicular to v23 in the plane defined by n1 and n2
    Coords3D m = v23.cross(n1).normalize();

    // Torsion angle calculation
    double x = n1.dot(n2);
    double y = m.dot(n2);
    double torsionAngle = atan2(y, x);
    cout << "torsion ange from the new method" << torsionAngle << endl;


    return t1234 ;
}



// Calculate torsion angle between four 3-d cartesian coordinates, we need torsion angle in rad to later on 
//calculate periodic torsion force and energy
pair<double, vector<Coords3D>> PeriodicTorsionForce::calculateForces(vector<Coords3D>& AP, const TorsionParameters& TP) {
    Coords3D u21 = (AP[TP.p2] - AP[TP.p1]) / atom_distance(AP[TP.p2] , AP[TP.p1]);
    Coords3D u23 = (AP[TP.p2] - AP[TP.p3]) / atom_distance(AP[TP.p2] , AP[TP.p3]);
    Coords3D u32 = (AP[TP.p3] - AP[TP.p2]) / atom_distance(AP[TP.p3] , AP[TP.p2]);
    Coords3D u34 = (AP[TP.p3] - AP[TP.p4]) / atom_distance(AP[TP.p3] , AP[TP.p4]);

    Coords3D n1 = plane_norm_vec(u21, u23);// n1 normal vector of plane 1
    Coords3D n2 = plane_norm_vec(u32, u34);// n2 normal vector of plane 2

    double dp = n1.dot(n2);//dp=dot product of two planes normal vectors
    dp = max(min(dp, 1.0), -1.0);//making sure dp which is the cos of torsion angle is between -1 and +1
    double sign = (n2.dot(u34) >= 0) ? 1.0 : -1.0;
    double torsionangle = sign * acos(dp); // t1234 is now in radians

    // energy derivative with respect to the torsion angle
    double dEdPhi = TP.k * TP.periodicity * sin(TP.periodicity * torsionangle - TP.phase);

    // calculate forces for each atom
    // the force is distributed among the atoms involved in the torsion according to their contribution to the torsion angle
    Coords3D v12 = u21;
    Coords3D v23 = u32;
    Coords3D v34 = -u34;
    double v23lengthsquared = v23.dot(v23);
    double LA1 = v12.cross(v23).dot(v12.cross(v23));//squared magnitude of the "lever arm" vector perpendicular to plane 1
    Coords3D a1 = n1 * (dEdPhi / LA1);
    double LA2 = v34.cross(v23).dot(v34.cross(v23));//squared magnitude of the "lever arm" vector perpendicular to plane 2
    Coords3D a4 = n2 * (dEdPhi / LA2);
    double v12projv23 = (-v12.dot(v23) / v23lengthsquared);//fraction of the v12 vector projection onto v23, normalized by the length of v23. The negative sign is used to adjust the direction of the force, ensuring it acts opposite to the projection of v12 onto v23.
    Coords3D a2 = (v23 * v12projv23).cross(a1);
    double v34projv23 = (-v34.dot(v23) / v23lengthsquared);//fraction of the v34 vector projection onto v23, normalized by the length of v23. The negative sign is used to adjust the direction of the force, ensuring it acts opposite to the projection of v34 onto v23.
    Coords3D a3 = (v23 * v34projv23).cross(a4);

    // forces on atoms are the negative of these, since force is the negative gradient of potential energy
    Coords3D f1 = -a1;//f1.x, f1.y, and f1.z are the force components acting on atom 1 in the x, y, and z directions, respectively.
    Coords3D f2 = -(a2 - a1);
    Coords3D f3 = -(a3 - a4);
    Coords3D f4 = -a4;

    vector<Coords3D> forces = { f1, f2, f3, f4 };
    return { torsionangle, forces };
}


//struct TorsionParameters {
//    int p1, p2, p3, p4;
//    double k, phase;
//    int periodicity;
//};

//Coords3D PeriodicTorsionForce::calculateForces(const Coords3D& atomPositions, const TorsionParameters& torsionParams, double torsionAngle) {
//    vector<vector<Coords3D>> forces(atomPositions.size(), vector<Coords3D>(1, Coords3D(0, 0, 0))); // Initialize forces to zero
//
//    // Simplified example to demonstrate concept - detailed force computation would be more complex
//    double deltaAngle = torsionParams.periodicity * torsionAngle - torsionParams.phase;
//    double dEdAngle = -torsionParams.k * torsionParams.periodicity * sin(deltaAngle);
//
//    //now we can assume that this torsion is devided 
//
//    // Placeholder for force vector calculation - in reality, you'd compute these based on geometry of the torsion
//    // For this example, assume each atom equally shares the force derived from dEdAngle, equally distributed among XYZ
//    Coords3D force(dEdAngle / 4.0, dEdAngle / 4.0, dEdAngle / 4.0); // Dividing by 4 as a placeholder for equal distribution
//
//    // Assign calculated force to each atom involved in the torsion
//    for (size_t i = 0; i < atomPositions.size(); ++i) {
//        forces[i] = force; // In a real calculation, forces on each atom would differ
//    }
//
//    return forces;
//}
//



//static vector<Coords3D> calculateForces(const Coords3D& r1, const Coords3D& r2, const Coords3D& r3, const Coords3D& r4, double k, double phase, int n) {
//    // Calculate vectors between atoms
//    Coords3D v12 = r2 - r1;
//    Coords3D v23 = r3 - r2;
//    Coords3D v34 = r4 - r3;
//
//    // Calculate normals to the planes defined by first three and last three atoms
//    Coords3D n1 = v12.cross(v23).normalize();
//    Coords3D n2 = v23.cross(v34).normalize();
//
//    // Calculate vector perpendicular to v23 in the plane defined by n1 and n2
//    Coords3D m = v23.cross(n1).normalize();
//
//    // Torsion angle calculation
//    double x = n1.dot(n2);
//    double y = m.dot(n2);
//    double torsionAngle = atan2(y, x);
//    cout << "torsion ange from the new method" << torsionAngle << endl;
//
//    // Energy derivative with respect to the torsion angle
//    double dEdPhi = k * n * sin(n * torsionAngle - phase);
//
//    // Calculate forces for each atom
//    // The force is distributed among the atoms involved in the torsion according to their contribution to the torsion angle
//    double v23LengthSquared = v23.dot(v23);
//    Coords3D a1 = n1 * (dEdPhi / (v12.cross(v23).dot(v12.cross(v23))));
//    Coords3D a4 = n2 * (dEdPhi / (v34.cross(v23).dot(v34.cross(v23))));
//
//    Coords3D a2 = (v23 * (-v12.dot(v23) / v23LengthSquared)).cross(a1);
//    Coords3D a3 = (v23 * (-v34.dot(v23) / v23LengthSquared)).cross(a4);
//
//    // Forces on atoms are the negative of these, since force is the negative gradient of potential energy
//    Coords3D f1 = -a1;//f1.x, f1.y, and f1.z are the force components acting on atom 1 in the x, y, and z directions, respectively.
//    Coords3D f2 = -(a2 - a1);
//    Coords3D f3 = -(a3 - a4);
//    Coords3D f4 = -a4;
//
//    vector<Coords3D> forces = { f1, f2, f3, f4 };
//    return forces;
//}



//int main() {
//    // Example usage with Coords3D
//    Coords3D coords1(1.0, 0.6, 0);
//    Coords3D coords2(0, 0, 0);
//    Coords3D coords3(1.57, 0, 0);
//    Coords3D coords4(2.07, 2.8, 0);
//    double torsionAngle = PeriodicTorsionForce::torsion_angle(coords1, coords2, coords3, coords4);
//    cout << "Torsion angle: " << torsionAngle << " degrees" << endl;
//    return 0;
//}


