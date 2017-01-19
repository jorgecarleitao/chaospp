#ifndef chaospp_auxiliar_h
#define chaospp_auxiliar_h

#include <utility> // for std::pair

#include <mpreal.h>
#include <Eigen/Dense>

typedef mpfr::mpreal Float;
using Matrix = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<Float, Eigen::Dynamic, 1>;


namespace aux {
    const Float pi = atan(Float(1))*4;

    typedef std::pair<Float, Float> pair;

    inline Float urandom() {
        return mpfr::random();
    }

    inline Float nrandom() {
        return mpfr::grandom();
    }

    Vector unitaryVector(unsigned int D) {
        Vector vector(D);
        Float norm = 0;
        for (unsigned int d = 0; d < D; d++) {
            vector[d] = nrandom();
            norm += vector[d]*vector[d];
        }
        norm = sqrt(norm);
        for (size_t d = 0; d < D; d++) {
            vector[d] /= norm;
        }
        return vector;
    }

    inline Float get_norm(Vector const& vector)
    {
        Float norm = 0;
        for(unsigned int i = 0; i < vector.size(); i++)
            norm += vector[i]*vector[i];
        return sqrt(norm);
    }

    inline void normalize(Vector & vector, Float const& norm)
    {
        for(unsigned int i = 0; i < vector.size(); i++)
            vector[i] /= norm;
    }

    inline void normalize(Vector & vector)
    {
        Float norm = get_norm(vector);
        normalize(vector, norm);
    }
}


#endif
