#include "NonbondedForce.h"
#include <iostream>  // Include this at the top of your file
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include "SystemXMLParser.h"



using namespace std;
using namespace BaseLine;

const double TWO_OVER_SQRT_PI = 1.1283791670955126; // Accurate value of 2/sqrt(pi)
//#define ONE_4PI_EPS0 8.9875517873681764e9 //​ unit N.m2/c2
#define ONE_4PI_EPS0 138.93545764446428693393 //kJ⋅nm/(mol⋅e2)
                     

NonbondedForce::NonbondedForce() {}

// Initialize reciprocal box vectors
void NonbondedForce::initializeReciprocalVectors(const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    //recipBoxVec[0] = 2 * PI / boxInfo.boxSize[0];
    //recipBoxVec[1] = 2 * PI / boxInfo.boxSize[1];
    //recipBoxVec[2] = 2 * PI / boxInfo.boxSize[2];
    for (int d = 0; d < 3; d++) {

        recipBoxVec[d] = 1 / boxInfo.boxSize[d];

        //recipBoxVec[0] = 1 / boxInfo.boxSize[0];
        //recipBoxVec[1] = 1 / boxInfo.boxSize[1];
        //recipBoxVec[2] = 1 / boxInfo.boxSize[2];
    }



}

// finding the first power of two above the box size which is in angstrom unit, roughly each A needs one grid point 
int NonbondedForce::findFFTDimension(int minimum) {
    // Calculate the number of grid points as one per angstrom
    //bug found: since the input from OpenMM is in nm and roughly each A needs one grid point, minimum should multiplied to 10 to found the number of grids
    minimum = minimum * 10;
    int power_of_two = 1;
    while (power_of_two < minimum) {
        power_of_two *= 2;
    }
    // Find the next highest power of two
    return power_of_two;
}

void NonbondedForce::initializeGridSize(const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    gridSize[0] = findFFTDimension(static_cast<int>(ceil(boxInfo.boxSize[0])));
    gridSize[1] = findFFTDimension(static_cast<int>(ceil(boxInfo.boxSize[1])));
    gridSize[2] = findFFTDimension(static_cast<int>(ceil(boxInfo.boxSize[2])));
}


void NonbondedForce::InitparticleIndexAndFraction(const vector<Coords3D>& atomPositions) {
    // Initialize particleIndex and particleFraction
    int numAtoms = atomPositions.size();
    particleIndex.resize(numAtoms, vector<int>(3, 0));
    particleFraction.resize(numAtoms, vector<double>(3, 0.0));

    // Compute particleIndex and particleFraction
    for (int i = 0; i < numAtoms; i++) {
        Coords3D coord = atomPositions[i];
        for (int d = 0; d < 3; d++) {
            double t = coord[d] * recipBoxVec[d];
            t = (t - floor(t)) * gridSize[d];
            int ti = static_cast<int>(t);
            particleFraction[i][d] = t - ti;
            particleIndex[i][d] = ti % gridSize[d];
        }
    }



}



void NonbondedForce::computeBSplineParameters(int pmeOrder, const vector<Coords3D>& atomPositions) {
    // Initialize the data array for B-spline values
    vector<double> data(pmeOrder, 0.0);
    vector<double> ddata(pmeOrder, 0.0);
    int maxSize = *max_element(gridSize.begin(), gridSize.end());
    vector<double> bsplines_data(maxSize, 0.0);
    int numAtoms = atomPositions.size();

    // Set initial conditions for B-splines
    data[pmeOrder - 1] = 0.0;
    data[1] = 0.0;
    data[0] = 1.0;

    // Compute B-spline values using recursive formula
    for (int i = 3; i < pmeOrder; i++) {
        double div = 1.0 / (i - 1);
        data[i - 1] = 0.0;
        for (int j = 1; j < (i - 1); j++) {
            data[i - j - 1] = div * (j * data[i - j - 2] + (i - j) * data[i - j - 1]);
        }
        data[0] = div * data[0];
    }

    // Differentiate B-spline values
    ddata[0] = -data[0];
    for (int i = 1; i < pmeOrder; i++) {
        ddata[i] = data[i - 1] - data[i];
    }

    // Normalize the last value
    double div = 1.0 / (pmeOrder - 1);
    data[pmeOrder - 1] = div * data[pmeOrder - 2];
    for (int i = 1; i < (pmeOrder - 1); i++) {
        data[pmeOrder - i - 1] = div * (i * data[pmeOrder - i - 2] + (pmeOrder - i) * data[pmeOrder - i - 1]);
    }
    data[0] = div * data[0];

    // Fill bsplines_data for Fourier transform
    for (int i = 1; i <= pmeOrder; i++) {
        bsplines_data[i] = data[i - 1];
    }

    // Resize the output vectors to match the grid size
    bsplines_moduli.resize(3);
    for (int d = 0; d < 3; d++) {
        bsplines_moduli[d].resize(gridSize[d], 0.0);
    }

    // Fourier transform to compute the moduli
    for (int d = 0; d < 3; d++) {
        for (int i = 0; i < gridSize[d]; i++) {
            double sc = 0.0;
            double ss = 0.0;
            for (int j = 0; j < gridSize[d]; j++) {
                double arg = (2 * PI * i * j) / gridSize[d];
                sc += bsplines_data[j] * cos(arg);
                ss += bsplines_data[j] * sin(arg);
            }
            bsplines_moduli[d][i] = sc * sc + ss * ss;
        }

        // Ensure no small or negative values in moduli
        for (int i = 0; i < gridSize[d]; i++) {
            if (bsplines_moduli[d][i] < 1e-7) {
                bsplines_moduli[d][i] = (bsplines_moduli[d][(i - 1 + gridSize[d]) % gridSize[d]] + bsplines_moduli[d][(i + 1) % gridSize[d]]) / 2.0;
            }
        }
    }

    // Resize bsplines_theta and bsplines_dtheta to match numAtoms and pmeOrder
    bsplines_theta.resize(3, vector<vector<double>>(numAtoms, vector<double>(pmeOrder, 0.0)));
    bsplines_dtheta.resize(3, vector<vector<double>>(numAtoms, vector<double>(pmeOrder, 0.0)));


    // Compute bsplines_theta and bsplines_dtheta for each particle and each dimension
    for (int i = 0; i < numAtoms; i++) {
        for (int d = 0; d < 3; d++) {
            double dr = particleFraction[i][d];

            bsplines_theta[d][i][pmeOrder - 1] = 0;
            bsplines_theta[d][i][1] = dr;
            bsplines_theta[d][i][0] = 1 - dr;

            for (int k = 3; k < pmeOrder; k++) {
                div = 1.0 / (k - 1.0);
                bsplines_theta[d][i][k - 1] = div * dr * bsplines_theta[d][i][k - 2];
                for (int l = 1; l < (k - 1); l++) {
                    bsplines_theta[d][i][k - l - 1] = div * ((dr + l) * bsplines_theta[d][i][k - l - 2] + (k - l - dr) * bsplines_theta[d][i][k - l - 1]);
                }
                bsplines_theta[d][i][0] = div * (1 - dr) * bsplines_theta[d][i][0];
            }

            // Differentiate
            bsplines_dtheta[d][i][0] = -bsplines_theta[d][i][0];

            for (int k = 1; k < pmeOrder; k++) {
                bsplines_dtheta[d][i][k] = bsplines_theta[d][i][k - 1] - bsplines_theta[d][i][k];
            }

            div = 1.0 / (pmeOrder - 1);
            bsplines_theta[d][i][pmeOrder - 1] = div * dr * bsplines_theta[d][i][pmeOrder - 2];

            for (int l = 1; l < (pmeOrder - 1); l++) {
                bsplines_theta[d][i][pmeOrder - l - 1] = div * ((dr + l) * bsplines_theta[d][i][pmeOrder - l - 2] + (pmeOrder - l - dr) * bsplines_theta[d][i][pmeOrder - l - 1]);
            }
            bsplines_theta[d][i][0] = div * (1 - dr) * bsplines_theta[d][i][0];
        }
    }
}




//void NonbondedForce::initializeBSplineModuli(int pmeOrder) {
//    // Compute moduli for X, Y, Z dimensions (assuming same gridSize for simplicity)
//    bsplines_moduli[0] = computeBSplineModuli(pmeOrder);
//    bsplines_moduli[1] = computeBSplineModuli(pmeOrder);
//    bsplines_moduli[2] = computeBSplineModuli(pmeOrder);
//}




void NonbondedForce::applyPeriodicBoundary(Coords3D& pos, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    pos[0] -= round(pos[0] / boxInfo.boxSize[0]) * boxInfo.boxSize[0];
    pos[1] -= round(pos[1] / boxInfo.boxSize[1]) * boxInfo.boxSize[1];
    pos[2] -= round(pos[2] / boxInfo.boxSize[2]) * boxInfo.boxSize[2];
}



void NonbondedForce::findAtomGridIndex(const vector<Coords3D>& atomPositions, vector<Int3D>& gridIndices, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];  // Use the half-complex format size

    for (size_t atom = 0; atom < atomPositions.size(); ++atom) {
        Coords3D pos = atomPositions[atom];
        applyPeriodicBoundary(pos, boxInfo);

        // Calculate the transformed position in reciprocal space and map it to the grid size
        double transformedX = (pos[0] * recipBoxVec[0] - floor(pos[0] * recipBoxVec[0])) * gridSizeX;
        double transformedY = (pos[1] * recipBoxVec[1] - floor(pos[1] * recipBoxVec[1])) * gridSizeY;
        double transformedZ = (pos[2] * recipBoxVec[2] - floor(pos[2] * recipBoxVec[2])) * gridSizeZ;

        // Create the grid index structure and push it back into the gridIndices vector
        Int3D index;
        index.x = static_cast<int>(transformedX) % gridSizeX;
        index.y = static_cast<int>(transformedY) % gridSizeY;
        index.z = static_cast<int>(transformedZ) % gridSizeZ;
        gridIndices.push_back(index);
    }
}



void NonbondedForce::gridSpreadCharge_Reference(const vector<Coords3D>& atomPositions, ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const int& pmeOrder) {
    int numAtoms = atomPositions.size();
    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];

    // Reset the grid
    fill(pmeGrid.begin(), pmeGrid.end(), Complex(0, 0));

    for (int i = 0; i < numAtoms; ++i) {
        double q = params.particles[i].q;

        // Grid index for the actual atom position
        int x0index = particleIndex[i][0];
        int y0index = particleIndex[i][1];
        int z0index = particleIndex[i][2];

        // Bspline factors for this atom in each dimension, calculated from fractional coordinates
        double* thetax = &bsplines_theta[0][i][0];
        double* thetay = &bsplines_theta[1][i][0];
        double* thetaz = &bsplines_theta[2][i][0];

        // Loop over norder*norder*norder (typically 5*5*5) neighbor cells
        for (int ix = 0; ix < pmeOrder; ++ix) {
            // Calculate index, apply PBC so we spread to index 0/1/2 when a particle is close to the upper limit of the grid
            int xindex = (x0index + ix) % gridSizeX;

            for (int iy = 0; iy < pmeOrder; ++iy) {
                int yindex = (y0index + iy) % gridSizeY;

                for (int iz = 0; iz < pmeOrder; ++iz) {
                    int zindex = (z0index + iz) % gridSizeZ;
                    // Calculate index in the charge grid
                    int index = xindex * gridSizeY * gridSizeZ + yindex * gridSizeZ + zindex;
                    // Add the charge times the bspline spread/interpolation factors to this grid position
                    pmeGrid[index] += q * thetax[ix] * thetay[iy] * thetaz[iz];
                    //cout << pmeGrid[index] << endl;
                }
            }
        }
    }
}


//void NonbondedForce::perform3DFFT(ComplexGrid& pmeGrid, bool forward) {
//    vector<size_t> shape = { static_cast<size_t>(gridSize[0]), static_cast<size_t>(gridSize[1]), static_cast<size_t>(gridSize[2]) };
//    vector<size_t> axes = { 0, 1, 2 };
//    vector<ptrdiff_t> stride = {
//        static_cast<ptrdiff_t>(gridSize[1] * gridSize[2] * sizeof(Complex)),
//        static_cast<ptrdiff_t>(gridSize[2] * sizeof(Complex)),
//        static_cast<ptrdiff_t>(sizeof(Complex))
//    };
//
//    pocketfft::c2c(shape, stride, stride, axes, forward, pmeGrid, pmeGrid, 1.0, 0);
//}

void NonbondedForce::perform3DFFT(ComplexGrid& pmeGrid, bool forward) {
    // Get the dimensions of the grid
    int nx = gridSize[0];
    int ny = gridSize[1];
    int nz = gridSize[2];

    // Create FFTW plan
    fftw_plan plan;
    if (forward) {
        plan = fftw_plan_dft_3d(nx, ny, nz,
            reinterpret_cast<fftw_complex*>(pmeGrid.data()),
            reinterpret_cast<fftw_complex*>(pmeGrid.data()),
            FFTW_FORWARD, FFTW_ESTIMATE);
    }
    else {
        plan = fftw_plan_dft_3d(nx, ny, nz,
            reinterpret_cast<fftw_complex*>(pmeGrid.data()),
            reinterpret_cast<fftw_complex*>(pmeGrid.data()),
            FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    // Execute the plan
    fftw_execute(plan);

    // Destroy the plan
    fftw_destroy_plan(plan);
}


void NonbondedForce::InitializeAlphaEwald(const NonbondedParams& params) {
    double tolerance = params.ewaldTolerance;
    double cutoff = params.cutoff;
    //double erfc_value;

    string methodType = "approximation";// the other option is "loop"

    if (methodType == "approximation") {
        alphaEwald = (1.0 / cutoff) * sqrt(-log(2.0 * tolerance));
    }
    else if(methodType == "loop") {
        double alpha_lo = 0.0;
        double alpha_hi = 10.0;
        double alpha = 0.0;

        for (int i = 0; i < 50; ++i) {
            alpha = (alpha_lo + alpha_hi) / 2.0;
            //erfc_value= erfc(alpha * cutoff);
            //cout << "Calculated erfc_value for alpha " <<alpha <<"is "<< erfc_value << "" << endl;
            if (erfc(alpha * cutoff) / cutoff >= tolerance) {
                alpha_lo = alpha;
            }
            else {
                alpha_hi = alpha;
            }
        }
        alphaEwald = alpha;//since cutoff in NexaBind is in nm unit alpha is a number between 0 to 10

    }




    // Debug output to verify the calculated value
    //cout << "Calculated alphaEwald: " << alphaEwald << endl;
}

void NonbondedForce::reciprocalConvolution_reference(ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& energy) {
    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];
    double reciprocalEnergy = 0.0;
    //int halfGridSizeZ = gridSizeZ / 2 + 1;


    //epsilon_r   Dielectric coefficient, typically 1.0.

    double one_4pi_eps = ONE_4PI_EPS0 / epsilon_r;
    double factor = PI * PI / (alphaEwald * alphaEwald);
    double boxfactor = PI * boxInfo.boxSize[0] * boxInfo.boxSize[1] * boxInfo.boxSize[2];


    //double boxfactor = -2 * PI * sqrt(PI) / (6.0 * boxInfo.boxSize[0] * boxInfo.boxSize[1] * boxInfo.boxSize[2]);
    double esum = 0.0;

    int maxkx = (gridSizeX + 1) / 2;
    int maxky = (gridSizeY + 1) / 2;
    int maxkz = (gridSizeZ + 1) / 2;

    //// Print the first 20 members of pmeGrid before the convolution
    //cout << "Before Convolution:" << endl;
    //for (int i = 0; i < 20 && i < pmeGrid.size(); ++i) {
    //    cout << "pmeGrid[" << i << "] = (" << pmeGrid[i].real() << ", " << pmeGrid[i].imag() << ")" << endl;
    //}

    double bfac = PI / alphaEwald;
    double fac1 = 2.0 * PI * PI * PI * sqrt(PI);
    double fac2 = alphaEwald * alphaEwald * alphaEwald;
    double fac3 = -2.0 * alphaEwald * PI * PI;

    for (int kx = 0; kx < gridSizeX; ++kx) {
        double mx = (kx < maxkx) ? kx : (kx - gridSizeX);
        double mhx = mx * recipBoxVec[0];
        double bx = boxfactor*bsplines_moduli[0][kx];

        for (int ky = 0; ky < gridSizeY; ++ky) {
            double my = (ky < maxky) ? ky : (ky - gridSizeY);
            double mhy = my * recipBoxVec[1];
            double by = bsplines_moduli[1][ky];

            for (int kz = 0; kz < gridSizeZ; ++kz) {
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }

                double mz = (kz < maxkz) ? kz : (kz - gridSizeZ);
                double mhz = mz * recipBoxVec[2];

                //int index = (kx * gridSizeY + ky) * halfGridSizeZ + kz;
                int index = kx * gridSizeY* gridSizeZ + ky * gridSizeZ + kz;

                Complex& gridValue = pmeGrid[index]; // Accessing the grid value

                double d1 = gridValue.real();
                double d2 = gridValue.imag();

                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                double bz = bsplines_moduli[2][kz];
                double denom = m2 * bx * by * bz;

                double eterm = one_4pi_eps * exp(-factor * m2) / denom;

                // Updating the pmeGrid value with the convolution term
                pmeGrid[index] = Complex(d1 * eterm, d2 * eterm);
                //cout << pmeGrid[index] << endl;

                double struct2 = d1 * d1 + d2 * d2;
                double ets2 = eterm * struct2;
                reciprocalEnergy += ets2;
            }
        }
    }

    energy = 0.5 * reciprocalEnergy; // Final energy calculation
}






void NonbondedForce::initializeErfcTable(const NonbondedParams& params) {
    // alphaEwald is already calculated
    float ewaldDX = params.cutoff / Num_Table_Points;  // Interval between points
    for (int i = 0; i < Num_Table_Points; ++i) {
        float r = i * ewaldDX;
        float alphaR = alphaEwald * r;
        erfcTable.push_back(erfc(alphaR));
    }
}



//new reference version
void NonbondedForce::ReciprocalForcesAndEnergy(vector<Coords3D>& forces, const vector<Coords3D>& atomPositions, const ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy, const int& pmeOrder) {
    int N = atomPositions.size();
    double directEnergy = 0.0;


    //// Direct space (short-range) part for forces and energy
    //for (int i = 0; i < N; ++i) {
    //    for (int j = i + 1; j < N; ++j) {
    //        Coords3D diff = atomPositions[j] - atomPositions[i];
    //        double r = sqrt(diff.dot(diff)); // Calculate the norm directly without normalizing
    //        Coords3D unitDiff = diff / r; // Normalize returns the unit vector, norm() gives its length

    //        if (r < params.cutoff) {
    //            double qiqj = params.particles[i].q * params.particles[j].q;
    //            double erfVal = erfcTable[static_cast<int>(r * erfcTable.size() / params.cutoff)];
    //            double forceMagnitude = qiqj * erfVal / (r * r);
    //            Coords3D forceVector = unitDiff * forceMagnitude; // normalize() to get the direction
    //            forces[i] += forceVector;
    //            forces[j] -= forceVector;
    //            directEnergy += qiqj * erfVal / r;
    //        }
    //    }
    //}

    // Reciprocal space (long-range) part for forces and energy from PME grid

    // Initialize variables
    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];

    // Interpolate forces from the PME grid
    for (int i = 0; i < atomPositions.size(); ++i) {
        double fx = 0.0, fy = 0.0, fz = 0.0;
        double q = params.particles[i].q;

        int x0index = particleIndex[i][0];
        int y0index = particleIndex[i][1];
        int z0index = particleIndex[i][2];

        for (int ix = 0; ix < pmeOrder; ix++) {
            int xindex = (x0index + ix) % gridSizeX;
            double tx = bsplines_theta[0][i][ix];
            double dtx = bsplines_dtheta[0][i][ix];

            for (int iy = 0; iy < pmeOrder; iy++) {
                int yindex = (y0index + iy) % gridSizeY;
                double ty = bsplines_theta[1][i][iy];
                double dty = bsplines_dtheta[1][i][iy];

                for (int iz = 0; iz < pmeOrder; iz++) {
                    int zindex = (z0index + iz) % gridSizeZ;
                    double tz = bsplines_theta[2][i][iz];
                    double dtz = bsplines_dtheta[2][i][iz];
                    int index = xindex * gridSizeY * gridSizeZ + yindex * gridSizeZ + zindex;

                    double gridvalue = pmeGrid[index].real();

                    fx += dtx * ty * tz * gridvalue;
                    fy += tx * dty * tz * gridvalue;
                    fz += tx * ty * dtz * gridvalue;
                }
            }
        }

        forces[i][0] -= q * (fx * gridSizeX * recipBoxVec[0]);
        forces[i][1] -= q * (fy * gridSizeY * recipBoxVec[1]);
        forces[i][2] -= q * (fz * gridSizeZ * recipBoxVec[2]);
    }


    // Combine both energies
    totalPEnergy += directEnergy;

    // Optional: Output results for debugging
    // cout << "Total Electrostatic Energy: " << totalPEnergy << endl;
    // for (int i = 0; i < N; ++i) {
    //     cout << "Force on particle " << i << ": (" << forces[i][0] << ", " << forces[i][1] << ", " << forces[i][2] << ")" << endl;
    // }
}



//void computeNeighborList(NeighborList& neighborList, const vector<Coords3D>& atomPositions, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double maxDistance, double minDistance, bool usePeriodic){
//
//
//    neighborList.pairs.clear();
//    int nAtoms = atomPositions.size();
//    double maxDistanceSquared = maxDistance * maxDistance;
//    double minDistanceSquared = minDistance * minDistance;
//
//    for (int i = 0; i < nAtoms - 1; ++i) {
//        for (int j = i + 1; j < nAtoms; ++j) {
//            Coords3D diff = PeriodicBoundaryCondition::minimumImageVector(atomPositions[i], atomPositions[j], boxInfo);
//            double distanceSquared=diff.dot(diff);
//            if (distanceSquared <= maxDistanceSquared && distanceSquared >= minDistanceSquared) {
//                neighborList.pairs.push_back(make_pair(i, j));
//            }
//        }
//    }
//}
//
//void NonbondedForce::DirectForcesAndEnergy(vector<Coords3D>& forces, const vector<Coords3D>& atomPositions, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy) {
//    int N = atomPositions.size();
//    double totalVdwEnergy = 0.0;
//    double totalRealSpaceEwaldEnergy = 0.0;
//    double switchValue = 1.0, switchDeriv = 0.0;
//
//    for (const auto& pair : neighborList) {
//        int ii = pair.first;
//        int jj = pair.second;
//
//        Coords3D deltaR = atomPositions[jj] - atomPositions[ii];
//        boxInfo.applyMinimumImage(deltaR);
//        double r2 = deltaR.dot(deltaR);
//        double r = sqrt(r2);
//        double inverseR = 1.0 / r;
//
//        // Switch function
//        if (useSwitch && r > params.switchingDistance) {
//            double t = (r - params.switchingDistance) / (params.cutoff - params.switchingDistance);
//            switchValue = 1 + t * t * t * (-10 + t * (15 - t * 6));
//            switchDeriv = t * t * (-30 + t * (60 - t * 30)) / (params.cutoff - params.switchingDistance);
//        }
//
//        // Electrostatic interaction
//        double alphaR = alphaEwald * r;
//        double dEdR = ONE_4PI_EPS0 * params.particles[ii].q * params.particles[jj].q * inverseR * inverseR * inverseR;
//        dEdR *= (erfc(alphaR) + 2 * alphaR * exp(-alphaR * alphaR) / SQRT_PI);
//
//        // van der Waals interaction
//        double sig = params.particles[ii].sig + params.particles[jj].sig;
//        double sig2 = sig * inverseR;
//        sig2 *= sig2;
//        double sig6 = sig2 * sig2 * sig2;
//        double eps = params.particles[ii].eps * params.particles[jj].eps;
//        dEdR += switchValue * eps * (12.0 * sig6 - 6.0) * sig6 * inverseR * inverseR;
//        double vdwEnergy = eps * (sig6 - 1.0) * sig6;
//
//        // Apply switching function to the van der Waals force
//        if (useSwitch) {
//            dEdR -= vdwEnergy * switchDeriv * inverseR;
//            vdwEnergy *= switchValue;
//        }
//
//        // Accumulate forces
//        Coords3D force = deltaR * dEdR;
//        forces[ii] += force;
//        forces[jj] -= force;
//
//        // Accumulate energies
//        totalRealSpaceEwaldEnergy += ONE_4PI_EPS0 * params.particles[ii].q * params.particles[jj].q * inverseR * erfc(alphaR);
//        totalVdwEnergy += vdwEnergy;
//    }
//
//    totalPEnergy += totalRealSpaceEwaldEnergy + totalVdwEnergy;
//}










vector<Coords3D> NonbondedForce::calculateForces(const vector<Coords3D>& atomPositions, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    // Initialize the PME grid as a complex grid
    
    initializeReciprocalVectors(boxInfo);
    initializeGridSize(boxInfo);
    const int pmeOrder = 5;//to have cubic b-spline
    //initializeBSplineModuli(pmeOrder);
    


    InitparticleIndexAndFraction(atomPositions);
    computeBSplineParameters(pmeOrder, atomPositions);



    ComplexGrid pmeGrid(gridSize[0] * gridSize[1] * gridSize[2], Complex(0.0, 0.0));
    // ComplexGrid pmeGrid(gridSize[0] * gridSize[1] * (gridSize[2] / 2 + 1), Complex(0.0, 0.0)); // bug resolved

    vector<Coords3D> forces(atomPositions.size(), { 0.0, 0.0, 0.0 });



    gridSpreadCharge_Reference(atomPositions, pmeGrid, params,boxInfo, pmeOrder);//here pmeGrid gets updated
    //reciprocalConvolution(pmeGrid, params);//here one more time pmeGrid gets updated

    // Perform 3D FFT
    perform3DFFT(pmeGrid, true);

    InitializeAlphaEwald(params);//required parameter for reciprocalConvolution_reference

    reciprocalConvolution_reference(pmeGrid, params,  boxInfo, totalPEnergy); // the reciprocal part of the potential energy is calculated here


    // Print the first 20 members of pmeGrid before the convolution
    cout << "After Convolution NexaBind:" << endl;
    for (int i = 0; i < 40 && i < pmeGrid.size(); ++i) {
        cout << "pmeGrid[" << i << "] = (" << pmeGrid[i].real() << ", " << pmeGrid[i].imag() << ")" << endl;
    }



    // Perform 3D FFT (inverse)
    perform3DFFT(pmeGrid, false);

    // Print the first 20 members of pmeGrid before the convolution
    cout << "After inverse NexaBind:" << endl;
    for (int i = 0; i < 40 && i < pmeGrid.size(); ++i) {
        cout << "pmeGrid[" << i << "] = (" << pmeGrid[i].real() << ", " << pmeGrid[i].imag() << ")" << endl;
    }

    // gridInterpolateForce(forces, atomPositions, pmeGrid, params, boxInfo);

    //update totalPEnergy
    //calculatePMEPotentialEnergy(pmeGrid, params, totalPEnergy);
    //DirectAndReciprocalForcesAndEnergy(forces, atomPositions,  pmeGrid, params, boxInfo, totalPEnergy, pmeOrder);


    initializeErfcTable(params);
    ReciprocalForcesAndEnergy(forces, atomPositions, pmeGrid, params, boxInfo, totalPEnergy, pmeOrder); // the direct part of the potential energy is calculated here and is added to the totalPEnergy

    // DirectForcesAndEnergy(forces, atomPositions, params, boxInfo, totalPEnergy);



    //applyExclusions(atomPositions, forces, totalPEnergy, params, boxInfo);

    return forces;
}