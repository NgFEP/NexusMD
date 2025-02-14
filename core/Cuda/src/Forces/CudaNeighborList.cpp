#include "CudaNeighborList.h"
#include <cmath>
#include <cassert>
#include <algorithm> // For min/max
#include <utility>   // For pair


using namespace std;
using namespace Cuda;


VoxelHash::VoxelHash(double vsx, double vsy, double vsz, const PeriodicBoundaryCondition::BoxInfo& boxInfo, bool usePeriodic)
    : _voxelSizeX(vsx), _voxelSizeY(vsy), _voxelSizeZ(vsz), _boxInfo(boxInfo), usePeriodic(usePeriodic) {
    if (usePeriodic) {
        nx = static_cast<int>(floor(_boxInfo.boxSize[0] / _voxelSizeX + 0.5));
        ny = static_cast<int>(floor(_boxInfo.boxSize[1] / _voxelSizeY + 0.5));
        nz = static_cast<int>(floor(_boxInfo.boxSize[2] / _voxelSizeZ + 0.5));
        _voxelSizeX = _boxInfo.boxSize[0] / nx;
        _voxelSizeY = _boxInfo.boxSize[1] / ny;
        _voxelSizeZ = _boxInfo.boxSize[2] / nz;
    }
}

VoxelIndex VoxelHash::getVoxelIndex(const Coords3D& location) const {
    Coords3D r = location;
    if (usePeriodic) {
        r[0] -= floor(r[0] / _boxInfo.boxSize[0]) * _boxInfo.boxSize[0];
        r[1] -= floor(r[1] / _boxInfo.boxSize[1]) * _boxInfo.boxSize[1];
        r[2] -= floor(r[2] / _boxInfo.boxSize[2]) * _boxInfo.boxSize[2];
    }
    int x = static_cast<int>(floor(r[0] / _voxelSizeX));
    int y = static_cast<int>(floor(r[1] / _voxelSizeY));
    int z = static_cast<int>(floor(r[2] / _voxelSizeZ));
    return VoxelIndex(x, y, z);
}

void VoxelHash::insert(int& item, const Coords3D& location) {
    VoxelIndex voxelIndex = getVoxelIndex(location);
    if (voxelMap.find(voxelIndex) == voxelMap.end())
        voxelMap[voxelIndex] = Voxel();
    Voxel& voxel = voxelMap[voxelIndex];
    voxel.push_back(make_pair(&location, item));
}

void VoxelHash::getNeighbors(
    NeighborList& neighbors,
    const pair<const Coords3D*, int>& referencePoint,
    const vector<set<int>>& exclusions,
    bool reportSymmetricPairs,
    double maxDistance,
    double minDistance) const {

    const int atomI = referencePoint.second; // atomI is the index of the reference atom.
    const Coords3D& locationI = *referencePoint.first;
    double maxDistanceSquared = maxDistance * maxDistance;
    double minDistanceSquared = minDistance * minDistance;

    int dIndexX = int(maxDistance / _voxelSizeX) + 1;// dIndexX, dIndexY, dIndexZ are the ranges in voxel indices to search for neighbors.
    int dIndexY = int(maxDistance / _voxelSizeY) + 1;
    int dIndexZ = int(maxDistance / _voxelSizeZ) + 1;
    VoxelIndex centerVoxelIndex = getVoxelIndex(locationI);
    int minz = centerVoxelIndex.z - dIndexZ;
    int maxz = centerVoxelIndex.z + dIndexZ;
    if (usePeriodic)
        maxz = min(maxz, minz + nz - 1);

    for (int z = minz; z <= maxz; ++z) {
        int boxz = static_cast<int>(floor(static_cast<float>(z) / nz));
        int miny = centerVoxelIndex.y - dIndexY;
        int maxy = centerVoxelIndex.y + dIndexY;
        if (usePeriodic) {
            // double yoffset = boxz * _boxInfo.boxSize[1] / _voxelSizeY;
            double yoffset = 0.0;
            miny -= static_cast<int>(ceil(yoffset));
            maxy -= static_cast<int>(floor(yoffset));
            maxy = min(maxy, miny + ny - 1);
        }
        for (int y = miny; y <= maxy; ++y) {
            int boxy = static_cast<int>(floor(static_cast<float>(y) / ny));
            int minx = centerVoxelIndex.x - dIndexX;
            int maxx = centerVoxelIndex.x + dIndexX;
            if (usePeriodic) {
                // double xoffset = (boxy * _boxInfo.boxSize[0] + boxz * _boxInfo.boxSize[0]) / _voxelSizeX;
                double xoffset = 0.0;
                minx -= static_cast<int>(ceil(xoffset));
                maxx -= static_cast<int>(floor(xoffset));
                maxx = min(maxx, minx + nx - 1);
            }
            for (int x = minx; x <= maxx; ++x) {
                VoxelIndex voxelIndex(x, y, z);
                if (usePeriodic) {
                    voxelIndex.x = (x + nx) % nx;
                    voxelIndex.y = (y + ny) % ny;
                    voxelIndex.z = (z + nz) % nz;
                }
                auto it = voxelMap.find(voxelIndex);// it as iterator
                if (it != voxelMap.end()) {
                    const Voxel& voxel = it->second;
                    for (const auto& entry : voxel) {
                        const Coords3D& locationJ = *entry.first;
                        int atomJ = entry.second;
                        if (atomI == atomJ)//It ignores the same atom (self-hits).
                            continue;
                        //Coords3D dr = locationJ - locationI;
                        //if (usePeriodic) {
                        //    dr[0] -= floor(dr[0] / _boxInfo.boxSize[0] + 0.5) * _boxInfo.boxSize[0];
                        //    dr[1] -= floor(dr[1] / _boxInfo.boxSize[1] + 0.5) * _boxInfo.boxSize[1];
                        //    dr[2] -= floor(dr[2] / _boxInfo.boxSize[2] + 0.5) * _boxInfo.boxSize[2];
                        //}
                        //double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

                        Coords3D dr = PeriodicBoundaryCondition::minimumImageVector(locationJ, locationI, _boxInfo);


                        double r2 = dr.dot(dr);

                        // Ignore exclusions.
                        if (exclusions[atomI].find(atomJ) != exclusions[atomI].end()) continue;


                        if (r2 >= minDistanceSquared && r2 < maxDistanceSquared) {
                            neighbors.addNeighbor(atomI, atomJ);
                            if (atomI < atomJ || reportSymmetricPairs) {
                                neighbors.addNeighbor(atomJ, atomI);
                                //if (exclusions.size() <= atomI || exclusions[atomI].find(atomJ) == exclusions[atomI].end())
                            }
                        }
                    }
                }
            }
        }
    }
}

void Cuda::obtainNeighborList(// Cuda:: should be mentioned for individual functions
    NeighborList& neighborList,
    const vector<Coords3D>& atomPositions,
    const PeriodicBoundaryCondition::BoxInfo& boxInfo,
    const vector<set<int>>& exclusions,
    double maxDistance,
    double minDistance,
    bool usePeriodic) {

    neighborList.clear();

    double edgeSizeX = maxDistance;
    double edgeSizeY = maxDistance;
    double edgeSizeZ = maxDistance;
    if (usePeriodic) {
        edgeSizeX = 0.5 * boxInfo.boxSize[0] / floor(boxInfo.boxSize[0] / maxDistance);
        edgeSizeY = 0.5 * boxInfo.boxSize[1] / floor(boxInfo.boxSize[1] / maxDistance);
        edgeSizeZ = 0.5 * boxInfo.boxSize[2] / floor(boxInfo.boxSize[2] / maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, boxInfo, usePeriodic); // Passing boxInfo directly
    for (int atomJ = 0; atomJ < atomPositions.size(); ++atomJ) {
        const Coords3D& location = atomPositions[atomJ];
        voxelHash.getNeighbors(neighborList, { &location, atomJ }, exclusions, false, maxDistance, minDistance);
        voxelHash.insert(atomJ, location);
    }
}
