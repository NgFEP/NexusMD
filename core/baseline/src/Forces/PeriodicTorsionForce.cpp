#include "stdafx.h"
#include "PeriodicTorsionForce.h"
#include "Coords3D.h"
#include <cmath>
#include <iomanip>//to print with 16 decimal  

#ifndef PI
#define PI 3.14159265358979323846
#endif

const double tenm3 = std::numeric_limits<double>::epsilon(); // Smallest difference between two doubles


using namespace std;
using namespace BaseLine; 

PeriodicTorsionForce::PeriodicTorsionForce() {
}

PeriodicTorsionForce::~PeriodicTorsionForce() {
}

// Calculate the distance between two points
double PeriodicTorsionForce::atom_distance(const Coords3D& coords1, const Coords3D& coords2) {
    Coords3D diff = coords2 - coords1;
    return sqrt(diff.dot(diff));
}


// Calculate torsion angle (in radians) between four points
double PeriodicTorsionForce::torsion_angle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const Coords3D& p4) {
    // Calculate vectors between atoms
    Coords3D v0 = p1 - p2;
    Coords3D v1 = p3 - p2;
    Coords3D v2 = p3 - p4;

    // Apply periodic boundary conditions if needed
    // APPLY_PERIODIC_TO_DELTA(v0)
    // APPLY_PERIODIC_TO_DELTA(v1)
    // APPLY_PERIODIC_TO_DELTA(v2)
    // Note: You need to define how to apply periodic conditions based on your simulation box

    // Cross products to get normals to planes formed by atoms
    Coords3D cp0 = v0.cross(v1);
    Coords3D cp1 = v1.cross(v2);


    //3 vectors between four atoms of i, j, k, l 
    //Coords3D rij = pi - pj;
    //Coords3D rkj = pk - pj;
    //Coords3D rkl = pk - pl;

    //Coords3D n1 = rij.cross(rkj);//normal vector of the ijk plane
    //Coords3D n2 = rkj.cross(rkl);//normal vector of the jkl plane

    //start of previous method
    Coords3D m1 = cp0.cross(cp1).normalize();

    double x = cp0.dot(cp1);
    double y = m1.dot(cp1);
    double angle2 = atan2(y, x);
    //end of previous method


    // Normalizing cross products to use for angle calculation
    double cosangle = cp0.normalize().dot(cp1.normalize());
    double angle;

    // Check if cosangle is near the singularity of acos()
    //if (cosangle > 0.99f || cosangle < -0.99f) {
    if (cosangle > 0.99f || cosangle < -0.99f) {// no f to make sure double precision is applied
        Coords3D cross_prod = cp0.cross(cp1);
        double scale = cp0.dot(cp0) * cp1.dot(cp1);
        angle = asin(sqrt(cross_prod.dot(cross_prod) / scale));
        if (cosangle < 0) {
            angle = PI - angle;
        }
    }
    else {
        angle = acos(cosangle);
    }

    // Adjust angle based on the direction relative to the v0 vector
    angle = (v0.dot(cp1) >= 0 ? angle : -angle);

    return angle;
}



vector<Coords3D> PeriodicTorsionForce::calculateForces(const vector<Coords3D>& atomPositions, const PTorsionParams& params) {
    // Assuming atomPositions contains positions for the four atoms involved in the torsion
    double torsionAngle = torsion_angle(atomPositions[0], atomPositions[1], atomPositions[2], atomPositions[3]);

    // Calculate energy derivative with respect to the torsion angle
    double deltaAngle = params.periodicity * torsionAngle - params.phase;
    // New line for energy calculation
    double energy = params.k * (1.0 + cos(deltaAngle));

    double dEdPhi = -params.k * params.periodicity * sin(deltaAngle);

    Coords3D r1 = atomPositions[0], r2 = atomPositions[1], r3 = atomPositions[2], r4 = atomPositions[3];

    //std::cout << std::fixed << std::setprecision(16) << "r1: " << r1 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "r2: " << r2 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "r3: " << r3 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "r4: " << r4 << std::endl;



    // Calculate vectors between atoms
    Coords3D v0 = r1 - r2;
    Coords3D v1 = r3 - r2;
    Coords3D v2 = r3 - r4;

    //std::cout << std::fixed << std::setprecision(20) << "v0: " << v0 << std::endl;
    //std::cout << std::fixed << std::setprecision(20) << "v1: " << v1 << std::endl;
    //std::cout << std::fixed << std::setprecision(20) << "v2: " << v2 << std::endl;


    // Calculate cross products
    Coords3D cp0 = v0.cross(v1);
    double uuu1= v0[1] * v1[2] - v0[2] * v1[1];
    double uuu2 = v0[1] * v1[2] ;
    //std::cout << std::fixed << std::setprecision(20) << "uuu2: " << uuu2 << std::endl;

    double uuu3 = - v0[2] * v1[1];

    Coords3D cp1 = v1.cross(v2);

    // Calculate norms of cross products and vectors
    double normCross1 = cp0.dot(cp0);
    double normSqrBC = v1.dot(v1);
    double normBC = sqrt(normSqrBC);
    double normCross2 = cp1.dot(cp1);
    double dp = 1.0 / normSqrBC;

    //for debugging purpose
    // Print norm values
    //std::cout << std::fixed << std::setprecision(16) << "cp0: (" << cp0[0] << ", " << cp0[1] << ", " << cp0[2] << ")" << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "cp1: (" << cp1[0] << ", " << cp1[1] << ", " << cp1[2] << ")" << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "theta: " << torsionAngle << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "dEdPhi: " << dEdPhi << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "normCross1: " << normCross1 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "normSqrBC: " << normSqrBC << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "normBC: " << normBC << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "normCross2: " << normCross2 << std::endl;
    ////std::cout << "dp: " << dp << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "dp: " << dp << std::endl;



    // Force factors calculation
    double factor1 = (-dEdPhi * normBC) / normCross1;
    double factor2 = v0.dot(v1) * dp;
    double factor3 = v2.dot(v1) * dp;
    double factor4 = (dEdPhi * normBC) / normCross2;

    // Print factor values
    //std::cout << std::fixed << std::setprecision(16) << "factor1: " << factor1 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "factor2: " << factor2 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "factor3: " << factor3 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "factor4: " << factor4 << std::endl;

    // Calculate forces based on cross products and factors
    Coords3D force1 = cp0 * factor1;
    Coords3D force4 = cp1 * factor4;
    Coords3D s = force1 * factor2 - force4 * factor3;
    //std::cout << std::fixed << std::setprecision(16) << "s: " << s << std::endl;

    Coords3D force2 = s - force1;
    Coords3D force3 = -s - force4;

    // Print factor values
    //std::cout << std::fixed << std::setprecision(16) << "force1: " << force1 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "force2: " << force2 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "force3: " << force3 << std::endl;
    //std::cout << std::fixed << std::setprecision(16) << "force4: " << force4 << std::endl;


     vector<Coords3D> forces = { force1, force2, force3, force4 };

    return forces;
}


//std::vector<Coords3D> PeriodicTorsionForce::calculateForces(const std::vector<Coords3D>& atomPositions, const TorsionParameters& params) {
//    double torsionAngle = torsion_angle(atomPositions[0], atomPositions[1], atomPositions[2], atomPositions[3]);
//    double dEdPhi = params.k * params.periodicity * sin(params.periodicity * torsionAngle - params.phase);
//
//    // Reusing vectors from the torsion_angle function for consistency
//    Coords3D v0 = atomPositions[0] - atomPositions[1];
//    Coords3D v1 = atomPositions[2] - atomPositions[1];
//    Coords3D v2 = atomPositions[2] - atomPositions[3];
//
//    Coords3D cp0 = v0.cross(v1);
//    Coords3D cp1 = v1.cross(v2);
//
//    double normCross1 = cp0.dot(cp0);
//    double normSqrBC = v1.dot(v1);
//    double normBC = std::sqrt(normSqrBC);
//    double normCross2 = cp1.dot(cp1);
//    double dp = 1.0 / normSqrBC;
//
//    Coords3D force1 = cp0 * ((-dEdPhi * normBC) / normCross1);
//    Coords3D force4 = cp1 * ((dEdPhi * normBC) / normCross2);
//    Coords3D s = force1 * (v0.dot(v1) * dp) - force4 * (v2.dot(v1) * dp);
//    Coords3D force2 = s - force1;
//    Coords3D force3 = -s - force4;
//
//    std::vector<Coords3D> forces = { force1, force2, force3, force4 };
//
//    return forces;
//}
