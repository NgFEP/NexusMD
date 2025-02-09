#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include "Coords3D.h"
#include "PeriodicBoundaryCondition.h"
#include <vector>
#include <set>
#include <map>
#include <utility>

namespace BaseLine {


    // To definie a type alias for a vector of pairs, representing a voxel
    using Voxel = std::vector<std::pair<const Coords3D*, int>>;


    class NeighborList {
    public:

        std::vector<std::pair<int, int>> pairs;

        void addNeighbor(int atomI, int atomJ) {
            pairs.push_back({ atomI, atomJ });
        }


        void clear() {
            pairs.clear();
        }

    };


    struct VoxelIndex {
        int x, y, z;
        VoxelIndex(int x, int y, int z) : x(x), y(y), z(z) {}
        bool operator<(const VoxelIndex& other) const {
            return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);
        }
    };

    class VoxelHash {
    public:
        VoxelHash(double vsx, double vsy, double vsz, const PeriodicBoundaryCondition::BoxInfo& boxInfo, bool usePeriodic);

        void insert(int& item, const Coords3D& location);
        VoxelIndex getVoxelIndex(const Coords3D& location) const;

        void getNeighbors(
            NeighborList& neighbors,
            const std::pair<const Coords3D*, int>& referencePoint,
            const std::vector<std::set<int>>& exclusions,
            bool reportSymmetricPairs,
            double maxDistance,
            double minDistance) const;

    private:
        double _voxelSizeX, _voxelSizeY, _voxelSizeZ;
        const PeriodicBoundaryCondition::BoxInfo _boxInfo; 
        int nx, ny, nz;
        const Coords3D* periodicBoxVectors;
        const bool usePeriodic;
        std::map<VoxelIndex, Voxel> voxelMap;
    };


    void obtainNeighborList(
        NeighborList& neighborList,
        const std::vector<Coords3D>& atomPositions,
        const PeriodicBoundaryCondition::BoxInfo& boxInfo,
        const std::vector<std::set<int>>& exclusions,
        double maxDistance,
        double minDistance,
        bool usePeriodic);

} // namespace BaseLine

#endif // NEIGHBORLIST_H
