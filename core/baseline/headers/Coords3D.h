#ifndef OPENMM_Coords3D_H_
#define OPENMM_Coords3D_H_


#include <cassert>
#include <iosfwd>
#include <cmath>

namespace BaseLine {

/**
 * This class represents a three component vector.  It is used for storing positions,
 * velocities, and forces.
 */

class Coords3D {
public:
    /**
     * Create a Coords3D whose elements are all 0.
     */
    Coords3D() {
        data[0] = data[1] = data[2] = 0.0;
    }
    /**
     * Create a Coords3D with specified x, y, and z components.
     */
    Coords3D(double x, double y, double z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }
    double& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    bool operator==(const Coords3D& rhs) const {
        return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
    }

    bool operator!=(const Coords3D& rhs) const {
        return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2]);
    }
    
    // Arithmetic operators
    
    // unary plus
    Coords3D operator+() const {
        return Coords3D(*this);
    }
    
    // plus
    Coords3D operator+(const Coords3D& rhs) const {
        const Coords3D& lhs = *this;
        return Coords3D(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }
    
    Coords3D& operator+=(const Coords3D& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    // unary minus
    Coords3D operator-() const {
        const Coords3D& lhs = *this;
        return Coords3D(-lhs[0], -lhs[1], -lhs[2]);
    }
    
    // minus
    Coords3D operator-(const Coords3D& rhs) const {
        const Coords3D& lhs = *this;
        return Coords3D(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
    }

    Coords3D& operator-=(const Coords3D& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

    // scalar product
    Coords3D operator*(double rhs) const {
        const Coords3D& lhs = *this;
        return Coords3D(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
    }

    Coords3D& operator*=(double rhs) {
        data[0] *= rhs;
        data[1] *= rhs;
        data[2] *= rhs;
        return *this;
    }

    // scalar division
    Coords3D operator/(double rhs) const {
        const Coords3D& lhs = *this;
        double scale = 1.0/rhs;
        return Coords3D(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale);
    }

    Coords3D& operator/=(double rhs) {
        double scale = 1.0/rhs;
        data[0] *= scale;
        data[1] *= scale;
        data[2] *= scale;
        return *this;
    }
    
    // dot product
    double dot(const Coords3D& rhs) const {
        const Coords3D& lhs = *this;
        return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
    }

    // cross product
    Coords3D cross(const Coords3D& rhs) const {
        return Coords3D(data[1]*rhs[2]-data[2]*rhs[1], data[2]*rhs[0]-data[0]*rhs[2], data[0]*rhs[1]-data[1]*rhs[0]);
    }


    // Existing constructor and other members...

    // Normalize the vector
    Coords3D normalize() {
        double len = std::sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
        Coords3D normalized(0.0, 0.0, 0.0); // Initialize a new Coords3D object
        if (len > 0) {
            normalized[0] = data[0] / len;
            normalized[1] = data[1] / len;
            normalized[2] = data[2] / len;
        }
        return normalized; // Return the new, normalized vector
    }

    
private:
    double data[3];
};

static Coords3D operator*(double lhs, Coords3D rhs) {
    return Coords3D(rhs[0]*lhs, rhs[1]*lhs, rhs[2]*lhs);
}

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Coords3D& v) {
    o<<'['<<v[0]<<", "<<v[1]<<", "<<v[2]<<']';
    return o;
}

} // namespace OpenMM

#endif /*OPENMM_Coords3D_H_*/
