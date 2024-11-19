#include "stdafx.h"
#include "CudaImproperTorsionForce.h"
#include <cmath>

using namespace std;
using namespace Cuda;


double ImproperTorsionForce::calculateImproperTorsionAngle(const Coords3D& p1, const Coords3D& p2, const Coords3D& p3, const Coords3D& p4) {
    // Calculate vectors between points
    Coords3D v1 = p3 - p1;
    Coords3D v2 = p3 - p2;
    Coords3D v3 = p4 - p3;

    // Calculate the bisector of v1 and v2
    Coords3D a = (v1 + v2) / 2.0;
    a = a.normalize();  // Normalize the bisector vector

    // Normalize vector v3
    Coords3D b = v3.normalize();

    // Calculate the dot product of vectors a (bisector) and b
    double cosangle = a.dot(b);

    double angle;

    // Check if cosangle is near the singularity of acos()
    if (cosangle > 0.99 || cosangle < -0.99) {
        Coords3D cross_prod = a.cross(b);
        double cross_mag = sqrt(cross_prod.dot(cross_prod));
        angle = asin(cross_mag);
        if (cosangle < 0) {
            angle = PI - angle;
        }
    }
    else {
        angle = acos(cosangle);
    }

    // Adjust angle based on the direction relative to the v1 vector
    // Coords3D cp1 = v1.cross(v2);
    // angle = (v1.dot(cp1) >= 0 ? angle : -angle);

    return angle; // Return angle in radians
}

vector<Coords3D> ImproperTorsionForce::calculateForces(const vector<Coords3D>& atomPositions, double k, double theta0) {
    double theta = calculateImproperTorsionAngle(atomPositions[0], atomPositions[1], atomPositions[2], atomPositions[3]);
    double deltaTheta = theta - theta0;
    double dEdTheta = k * deltaTheta;

    Coords3D v1 = atomPositions[1] - atomPositions[0];
    Coords3D v2 = atomPositions[2] - atomPositions[1];
    Coords3D v3 = atomPositions[3] - atomPositions[2];
    Coords3D n1 = v1.cross(v2);
    double rp = sqrt(n1.dot(n1));

    // Ensure the normal vector is non-zero
    if (rp == 0) return vector<Coords3D>(4, { 0, 0, 0 });

    Coords3D force_on_p4 = n1 * (dEdTheta / rp);

    // Applying equal and opposite forces to maintain the system's balance
    vector<Coords3D> forces(4);
    forces[0] = force_on_p4 * -0.25;
    forces[1] = force_on_p4 * -0.25;
    forces[2] = force_on_p4 * -0.25;
    forces[3] = force_on_p4 * 0.75; // Most of the corrective force applied to the atom most out of plane

    return forces;
}
