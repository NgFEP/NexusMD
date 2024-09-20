#ifndef CUDATypes_H_
#define CUDATypes_H_

#include <complex>
#include <vector>

namespace Cuda {

    struct Int3D
    {
        int x, y, z;
    };

    struct Double3D
    {
        double x, y, z;
    };

    // Define a complex grid type using std::complex for ease of use with FFTW
    typedef std::complex<double> Complex;
    typedef std::vector<Complex> ComplexGrid;


}
#endif /*CUDATypes_H_*/