#ifndef Types_H_
#define Types_H_

#include <complex>
#include <vector>

namespace BaseLine {

    struct Int3D
    {
        int x, y, z;
    };

    struct Double3D
    {
        double x, y, z;
    };

    // To define a complex grid type using std::complex for ease of use with FFTW
    typedef std::complex<double> Complex;
    typedef std::vector<Complex> ComplexGrid;


}
#endif /*Types_H_*/