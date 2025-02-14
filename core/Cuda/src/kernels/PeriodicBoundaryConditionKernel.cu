#include "PeriodicBoundaryConditionKernel.h" 

__device__ void minimumImageVectorDevice(double3* coords1, double3* coords2, double3* delta, double3* d_boxsize) {
    delta->x = coords2->x - coords1->x;
    delta->y = coords2->y - coords1->y;
    delta->z = coords2->z - coords1->z;

    if (delta->x > 0.5f * d_boxsize->x) {
        delta->x -= d_boxsize->x;
    }
    else if (delta->x < -0.5f * d_boxsize->x) {
        delta->x += d_boxsize->x;
    }

    if (delta->y > 0.5f * d_boxsize->y) {
        delta->y -= d_boxsize->y;
    }
    else if (delta->y < -0.5f * d_boxsize->y) {
        delta->y += d_boxsize->y;
    }

    if (delta->z > 0.5f * d_boxsize->z) {
        delta->z -= d_boxsize->z;
    }
    else if (delta->z < -0.5f * d_boxsize->z) {
        delta->z += d_boxsize->z;
    }
}
