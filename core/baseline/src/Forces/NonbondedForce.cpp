#include "NonbondedForce.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include "SystemXMLParser.h"



using namespace std;
using namespace BaseLine;



NonbondedForce::NonbondedForce() {}

// To initialize reciprocal box vectors
void NonbondedForce::initializeReciprocalVectors(const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    for (int d = 0; d < 3; d++) {

        recipBoxVec[d] = 1 / boxInfo.boxSize[d];


    }



}

// finding the first power of two above the box size which is in angstrom unit, roughly each A needs one grid point 
int NonbondedForce::findFFTDimension(int minimum) {
    // To calculate the number of grid points as one per angstrom
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
    // To initialize particleIndex and particleFraction
    int numAtoms = atomPositions.size();
    particleIndex.resize(numAtoms, vector<int>(3, 0));
    particleFraction.resize(numAtoms, vector<double>(3, 0.0));

    // To compute particleIndex and particleFraction
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
    // To initialize the data array for B-spline values
    vector<double> data(pmeOrder, 0.0);
    vector<double> ddata(pmeOrder, 0.0);
    int maxSize = *max_element(gridSize.begin(), gridSize.end());
    vector<double> bsplines_data(maxSize, 0.0);
    int numAtoms = atomPositions.size();

    // To set initial conditions for B-splines
    data[pmeOrder - 1] = 0.0;
    data[1] = 0.0;
    data[0] = 1.0;

    // To compute B-spline values using recursive formula
    for (int i = 3; i < pmeOrder; i++) {
        double div = 1.0 / (i - 1);
        data[i - 1] = 0.0;
        for (int j = 1; j < (i - 1); j++) {
            data[i - j - 1] = div * (j * data[i - j - 2] + (i - j) * data[i - j - 1]);
        }
        data[0] = div * data[0];
    }

    // To differentiate B-spline values
    ddata[0] = -data[0];
    for (int i = 1; i < pmeOrder; i++) {
        ddata[i] = data[i - 1] - data[i];
    }

    // To normalize the last value
    double div = 1.0 / (pmeOrder - 1);
    data[pmeOrder - 1] = div * data[pmeOrder - 2];
    for (int i = 1; i < (pmeOrder - 1); i++) {
        data[pmeOrder - i - 1] = div * (i * data[pmeOrder - i - 2] + (pmeOrder - i) * data[pmeOrder - i - 1]);
    }
    data[0] = div * data[0];

    // To fill bsplines_data for Fourier transform
    for (int i = 1; i <= pmeOrder; i++) {
        bsplines_data[i] = data[i - 1];
    }

    // To resize the output vectors to match the grid size
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

        // To ensure no small or negative values in moduli
        for (int i = 0; i < gridSize[d]; i++) {
            if (bsplines_moduli[d][i] < 1e-7) {
                bsplines_moduli[d][i] = (bsplines_moduli[d][(i - 1 + gridSize[d]) % gridSize[d]] + bsplines_moduli[d][(i + 1) % gridSize[d]]) / 2.0;
            }
        }
    }

    // To resize bsplines_theta and bsplines_dtheta to match numAtoms and pmeOrder
    bsplines_theta.resize(3, vector<vector<double>>(numAtoms, vector<double>(pmeOrder, 0.0)));
    bsplines_dtheta.resize(3, vector<vector<double>>(numAtoms, vector<double>(pmeOrder, 0.0)));


    // To compute bsplines_theta and bsplines_dtheta for each particle and each dimension
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

            // To differentiate
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


void NonbondedForce::applyPeriodicBoundary(Coords3D& pos, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {
    pos[0] -= round(pos[0] / boxInfo.boxSize[0]) * boxInfo.boxSize[0];
    pos[1] -= round(pos[1] / boxInfo.boxSize[1]) * boxInfo.boxSize[1];
    pos[2] -= round(pos[2] / boxInfo.boxSize[2]) * boxInfo.boxSize[2];
}



void NonbondedForce::findAtomGridIndex(const vector<Coords3D>& atomPositions, vector<Int3D>& gridIndices, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo) {

    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];

    for (int atom = 0; atom < atomPositions.size(); ++atom) {
        Coords3D pos = atomPositions[atom];
        applyPeriodicBoundary(pos, boxInfo);

        // To calculate the transformed position in reciprocal space and map it to the grid size
        double transformedX = (pos[0] * recipBoxVec[0] - floor(pos[0] * recipBoxVec[0])) * gridSizeX;
        double transformedY = (pos[1] * recipBoxVec[1] - floor(pos[1] * recipBoxVec[1])) * gridSizeY;
        double transformedZ = (pos[2] * recipBoxVec[2] - floor(pos[2] * recipBoxVec[2])) * gridSizeZ;

        // To create the grid index structure and push it back into the gridIndices vector
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

    // To reset the grid
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

        // To loop over norder*norder*norder (typically 5*5*5) neighbor cells
        for (int ix = 0; ix < pmeOrder; ++ix) {
            // To calculate index, apply PBC so we spread to index 0/1/2 when a particle is close to the upper limit of the grid
            int xindex = (x0index + ix) % gridSizeX;

            for (int iy = 0; iy < pmeOrder; ++iy) {
                int yindex = (y0index + iy) % gridSizeY;

                for (int iz = 0; iz < pmeOrder; ++iz) {
                    int zindex = (z0index + iz) % gridSizeZ;
                    // To calculate index in the charge grid
                    int index = xindex * gridSizeY * gridSizeZ + yindex * gridSizeZ + zindex;
                    // To add the charge times the bspline spread/interpolation factors to this grid position
                    pmeGrid[index] += q * thetax[ix] * thetay[iy] * thetaz[iz];
                }
            }
        }
    }
}



void NonbondedForce::perform3DFFT(ComplexGrid& pmeGrid, bool forward) {
    // To get the dimensions of the grid
    int nx = gridSize[0];
    int ny = gridSize[1];
    int nz = gridSize[2];

    // To create FFTW plan
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

    // To execute the plan
    fftw_execute(plan);

    // To destroy the plan
    fftw_destroy_plan(plan);
}


void NonbondedForce::InitializeAlphaEwald(const NonbondedParams& params) {
    double tolerance = params.ewaldTolerance;
    double cutoff = params.cutoff;

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

            if (erfc(alpha * cutoff) / cutoff >= tolerance) {
                alpha_lo = alpha;
            }
            else {
                alpha_hi = alpha;
            }
        }
        alphaEwald = alpha;//since cutoff in Nexus is in nm unit alpha is a number between 0 to 10

    }

}

void NonbondedForce::reciprocalConvolutionAndEnergy(ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& energy) {
    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];
    double reciprocalEnergy = 0.0;


    double one_4pi_eps = ONE_4PI_EPS0 / epsilon_r;
    double factor = PI * PI / (alphaEwald * alphaEwald);
    double boxfactor = PI * boxInfo.boxSize[0] * boxInfo.boxSize[1] * boxInfo.boxSize[2];

    double esum = 0.0;

    int maxkx = (gridSizeX + 1) / 2;
    int maxky = (gridSizeY + 1) / 2;
    int maxkz = (gridSizeZ + 1) / 2;


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

                int index = kx * gridSizeY* gridSizeZ + ky * gridSizeZ + kz;

                Complex& gridValue = pmeGrid[index]; 

                double d1 = gridValue.real();
                double d2 = gridValue.imag();

                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                double bz = bsplines_moduli[2][kz];
                double denom = m2 * bx * by * bz;

                double eterm = one_4pi_eps * exp(-factor * m2) / denom;

                // Updating the pmeGrid value with the convolution term
                pmeGrid[index] = Complex(d1 * eterm, d2 * eterm);

                double struct2 = d1 * d1 + d2 * d2;
                double ets2 = eterm * struct2;
                reciprocalEnergy += ets2;
            }
        }
    }

    energy = 0.5 * reciprocalEnergy;
}






void NonbondedForce::initializeErfcTable(const NonbondedParams& params) {
    // alphaEwald is already calculated
    double ewaldDX = params.cutoff / Num_Table_Points; 
    for (int i = 0; i < Num_Table_Points; ++i) {
        double r = i * ewaldDX;
        double alphaR = alphaEwald * r;
        erfcTable.push_back(erfc(alphaR));
    }
}


//new reference version
void NonbondedForce::PMEGridForces(vector<Coords3D>& forces, const vector<Coords3D>& atomPositions, const ComplexGrid& pmeGrid, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, double& totalPEnergy, const int& pmeOrder) {
    int N = atomPositions.size();
    double directEnergy = 0.0;

    // Reciprocal space (long-range) part for forces and energy from PME grid

    // To initialize variables
    int gridSizeX = gridSize[0];
    int gridSizeY = gridSize[1];
    int gridSizeZ = gridSize[2];

    // To interpolate forces from the PME grid
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


    // To combine both energies
    totalPEnergy += directEnergy;

}





void NonbondedForce::ReciprocalPMEcalculateForcesAndEnergy(vector<Coords3D>&forces, const vector<Coords3D>&atomPositions, double& totalPEnergy, const NonbondedParams & params, const PeriodicBoundaryCondition::BoxInfo & boxInfo) {

    
    // Initialize the PME grid as a complex grid
    
    initializeReciprocalVectors(boxInfo);
    initializeGridSize(boxInfo);
    const int pmeOrder = 5;// for cubic b-spline
    
    InitparticleIndexAndFraction(atomPositions);
    computeBSplineParameters(pmeOrder, atomPositions);

    ComplexGrid pmeGrid(gridSize[0] * gridSize[1] * gridSize[2], Complex(0.0, 0.0));


    gridSpreadCharge_Reference(atomPositions, pmeGrid, params,boxInfo, pmeOrder);//here pmeGrid gets updated

    // To perform 3D FFT
    perform3DFFT(pmeGrid, true);

    InitializeAlphaEwald(params);//required parameter for reciprocalConvolution_reference

    reciprocalConvolutionAndEnergy(pmeGrid, params,  boxInfo, totalPEnergy); // the reciprocal part of the potential energy is calculated here

    // Perform 3D FFT (inverse)
    perform3DFFT(pmeGrid, false);


    initializeErfcTable(params);
    PMEGridForces(forces, atomPositions, pmeGrid, params, boxInfo, totalPEnergy, pmeOrder); // the direct part of the potential energy is calculated here and is added to the totalPEnergy


}

// Direct Space Calculations


void NonbondedForce::NonbondedParamsModifier(const vector<Coords3D>& atomPositions, const NonbondedParams& params) {
    int numAtoms = atomPositions.size();
    _ModifiedNonbondedParams.particles.resize(numAtoms);
    
    for (int i = 0; i < numAtoms; i++) {
        _ModifiedNonbondedParams.particles[i].sig = 0.5 * params.particles[i].sig;
        _ModifiedNonbondedParams.particles[i].eps = 2.0 * sqrt(params.particles[i].eps);
        _ModifiedNonbondedParams.particles[i].q = params.particles[i].q;
    }
}

void NonbondedForce::CalculateSelfEnergy(const NonbondedParams& params, double& totalPEnergy) {
    double totalSelfEwaldEnergy = 0.0;

    int numberOfAtoms = params.particles.size();

    for (int atomID = 0; atomID < numberOfAtoms; atomID++) {
        double selfEwaldEnergy = ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[atomID].q * _ModifiedNonbondedParams.particles[atomID].q * alphaEwald / SQRT_PI;

        totalSelfEwaldEnergy -= selfEwaldEnergy;
    }
    totalPEnergy += totalSelfEwaldEnergy;

}


void NonbondedForce::DirectForcesAndEnergy(vector<Coords3D>& forces, const vector<Coords3D>& atomPositions, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const vector<set<int>>& exclusions, double& totalPEnergy) {
    NeighborList neighborList;
    double maxDistance = params.cutoff;
    double minDistance = 0.0;
    bool usePeriodic = true;

    obtainNeighborList(neighborList, atomPositions, boxInfo, exclusions, maxDistance, minDistance, usePeriodic);

    int N = atomPositions.size();
    double totalVdwEnergy = 0.0;
    double totalRealSpaceEwaldEnergy = 0.0;
    double switchValue = 1.0, switchDeriv = 0.0;

    for (const auto& pair : neighborList.pairs) {
        int ii = pair.first;
        //int ii = 7;
        int jj = pair.second;

        Coords3D deltaR = PeriodicBoundaryCondition::minimumImageVector(atomPositions[ii], atomPositions[jj], boxInfo);

        double r = deltaR.length();
        double inverseR = 1.0 / r;

        if (params.useSwitchingFunction && r > params.switchingDistance) {
            double t = (r - params.switchingDistance) / (params.cutoff - params.switchingDistance);
            switchValue = 1 + t * t * t * (-10 + t * (15 - t * 6));
            switchDeriv = t * t * (-30 + t * (60 - t * 30)) / (params.cutoff - params.switchingDistance);
        }

        double alphaR = alphaEwald * r;
        double dEdR = ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[ii].q * _ModifiedNonbondedParams.particles[jj].q * inverseR * inverseR * inverseR;
        dEdR *= (erfc(alphaR) + 2 * alphaR * exp(-alphaR * alphaR) / SQRT_PI);

        double sig = _ModifiedNonbondedParams.particles[ii].sig + _ModifiedNonbondedParams.particles[jj].sig;
        double sig2 = sig * inverseR;
        sig2 *= sig2;
        double sig6 = sig2 * sig2 * sig2;
        double eps = _ModifiedNonbondedParams.particles[ii].eps * _ModifiedNonbondedParams.particles[jj].eps;
        dEdR += switchValue * eps * (12.0 * sig6 - 6.0) * sig6 * inverseR * inverseR;
        double vdwEnergy = eps * (sig6 - 1.0) * sig6;

        if (params.useSwitchingFunction) {
            dEdR -= vdwEnergy * switchDeriv * inverseR;
            vdwEnergy *= switchValue;
        }

        Coords3D force = deltaR * dEdR;
        forces[ii] += force;
        forces[jj] -= force;

        totalRealSpaceEwaldEnergy += ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[ii].q * _ModifiedNonbondedParams.particles[jj].q * inverseR * erfc(alphaR);
        totalVdwEnergy += vdwEnergy;
    }

    totalPEnergy += totalRealSpaceEwaldEnergy + totalVdwEnergy;
}




void NonbondedForce::CalculateExclusionEnergyAndForces(vector<Coords3D>& forces, const vector<Coords3D>& atomPositions,  const PeriodicBoundaryCondition::BoxInfo& boxInfo, const vector<set<int>>& exclusions, double& totalPEnergy) {

    double totalExclusionEnergy = 0.0;

    int numberOfAtoms = atomPositions.size();

    for (int i = 0; i < numberOfAtoms; i++) {
        for (int exclusion : exclusions[i]) {
            if (exclusion > i) {
                int ii = i;
                int jj = exclusion;

                Coords3D deltaR = PeriodicBoundaryCondition::minimumImageVector(atomPositions[ii], atomPositions[jj], boxInfo);

                double r = deltaR.length(); 
                double inverseR = 1.0 / r;
                double alphaR = alphaEwald * r;
                double ExclusionEnergy = 0.0;

                if (erf(alphaR) > 1e-6) { // very small values close to zero can lead to numerical instability or insignificant contributions,therefore they are effectively ignored here
                    double dEdR = ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[ii].q * _ModifiedNonbondedParams.particles[jj].q * inverseR * inverseR * inverseR;
                    dEdR = dEdR * (erf(alphaR) - 2 * alphaR * exp(-alphaR * alphaR) / SQRT_PI);

                    Coords3D force = deltaR * dEdR;
                    forces[ii] -= force;
                    forces[jj] += force;

                    ExclusionEnergy = ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[ii].q * _ModifiedNonbondedParams.particles[jj].q * inverseR * erf(alphaR);
                }
                else {
                    ExclusionEnergy = alphaEwald * TWO_OVER_SQRT_PI * ONE_4PI_EPS0 * _ModifiedNonbondedParams.particles[ii].q * _ModifiedNonbondedParams.particles[jj].q;
                }


                totalExclusionEnergy += ExclusionEnergy;
            }
        }
    }
    totalPEnergy -= totalExclusionEnergy;

}


vector<Coords3D> NonbondedForce::calculateForces(const vector<Coords3D>& atomPositions, double& totalPEnergy, const NonbondedParams& params, const PeriodicBoundaryCondition::BoxInfo& boxInfo, const vector<set<int>>& exclusions) {
    vector<Coords3D> forces(atomPositions.size(), { 0.0, 0.0, 0.0 });

    ReciprocalPMEcalculateForcesAndEnergy(forces, atomPositions, totalPEnergy, params, boxInfo);

    NonbondedParamsModifier(atomPositions, params);

    CalculateSelfEnergy(params, totalPEnergy);

    DirectForcesAndEnergy(forces, atomPositions, params, boxInfo, exclusions, totalPEnergy);

    CalculateExclusionEnergyAndForces(forces, atomPositions, boxInfo, exclusions, totalPEnergy);

    return forces;
}