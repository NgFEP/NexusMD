#ifndef PERIODICTORSIONFORCE_H
#define PERIODICTORSIONFORCE_H

#include "Coords3D.h"
#include <vector>
#include "SystemXMLParser.h" //because of its Torsionparameters


namespace BaseLine {

    class PeriodicTorsionForce {
    public:
        PeriodicTorsionForce();
        ~PeriodicTorsionForce();

        // Calculate the distance between two 3D coordinates

        static double atom_distance(const Coords3D& coords1, const Coords3D& coords2);
        //static allows it to be called without an instance of the class.
        static Coords3D plane_norm_vec(const Coords3D& uvec1, const Coords3D& uvec2);

        // Calculate the torsion angle given four 3D coordinates
        static double torsion_angle(const Coords3D& coords1, const Coords3D& coords2, const Coords3D& coords3, const Coords3D& coords4);
        
        // Calculate the torsion angle and forces on all 4 involved atoms given four 3D coordinates
        static vector<Coords3D> calculateForces(vector<Coords3D>& AP, const TorsionParameters& TP);

    };

} // namespace BaseLine

#endif // PERIODICTORSIONFORCE_H