#ifndef __chaospp__map__
#define __chaospp__map__

#include <vector>
#include "auxiliar.h"
#include "io.h"


namespace map {


//! A general class of a map.
class Map {
public:
    const unsigned int D;
    std::string name;
    std::vector<aux::pair> boundary;
protected:
    Matrix _jacobian;  // stores the jacobian to avoid repeating allocations.

public:
    //! returns the point inside the boundary conditions
    static void apply_boundary_conditions(Vector &point, std::vector<aux::pair> const& bounding_box) {
        for(unsigned int i = 0; i < point.size(); i++) {
            while(point[i] > bounding_box[i].second)
                point[i] -= bounding_box[i].second - bounding_box[i].first;
            while(point[i] < bounding_box[i].first)
                point[i] += bounding_box[i].second - bounding_box[i].first;
        }
    }

    Map(unsigned int D, std::string name) : D(D), _jacobian(D,D), name(name), boundary(D) {}

    //! one time evolution of the map
    virtual void T(Vector & point) = 0;

    //! one time evolution of the map on the tangent space
    inline void dT(Vector const& point, Vector & vector) {
        vector = jacobian(point)*vector;
    }

    //! one time evolution of the map on the tangent space
    inline void dTmatrix(Vector const& point, Matrix & matrix) {
        matrix = jacobian(point)*matrix;
    }

    //! the jacobian matrix of the map at a given point.
    virtual Matrix const& jacobian(Vector const& point) = 0;

    //! test if point has left the restraining region.
    //! optional method in case you don't want to even use escape time
    virtual bool has_exited(Vector const& point) {return true;}
};


class Manneville : public Map {
    Float z;
    std::vector<aux::pair> bounding_box;
public:
    Manneville(Float z = 2) : Map(1, format("pm%1.f", z.toDouble())), bounding_box(1), z(z) {
        assert(z >= 1); // z = 1 is the Bernoulli shift.
        bounding_box[0] = aux::pair(0, 1);
        boundary[0] = aux::pair(0, 1);
    }

    void T(Vector & point) {
        point[0] = pow(point[0], z) + point[0];
        apply_boundary_conditions(point, bounding_box);
    }

    Matrix const& jacobian(Vector const& point) {
        _jacobian(0,0) = 1 + pow(point[0], z - 1)*z;
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (point[0] > Float("0.8"))
            return true;
        else
            return false;
    }
};


class Standard : public Map {
    Float k;
    std::vector<aux::pair> bounding_box;
public:
    Standard(Float k) : Map(2, format("sm%1.f", k.toDouble())), bounding_box(2), k(k/(2*aux::pi)) {
        bounding_box[0] = aux::pair(0, 1);
        bounding_box[1] = aux::pair(0, 1);

        boundary[0] = aux::pair(0, 1);
        boundary[1] = aux::pair(0, 1);
    }

    void T(Vector & point) {
        point[0] += k*sin(2*aux::pi*point[1]);
        point[1] += point[0];
        apply_boundary_conditions(point, bounding_box);
    }

    Matrix const& jacobian(Vector const& point) {
        _jacobian(0,0) = 1;
        _jacobian(1,0) = 1;
        _jacobian(0,1) = k*2*aux::pi*cos(2*aux::pi*point[1]);
        _jacobian(1,1) = 1 + _jacobian(0,1);
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (point[1] < 0.1)
            return true;
        else
            return false;
    };
};


//! Eq. (1) in http://arxiv.org/pdf/1311.7632v2.pdf
class CoupledStandard : public Map {
    Float k1, k2, xi;
    std::vector<aux::pair> bounding_box;
public:
    CoupledStandard() : Map(4, "csm"), bounding_box(4), k1("-2.25"/(2*aux::pi)), k2("-3.0"/(2*aux::pi)), xi("1.0"/(2*aux::pi)) {
        //! right after Eq. (1)
        bounding_box[0] = aux::pair("-0.5", "0.5");
        bounding_box[1] = aux::pair("-0.5", "0.5");
        bounding_box[2] = aux::pair("-0.5", "0.5");
        bounding_box[3] = aux::pair("-0.5", "0.5");

        for (unsigned int i = 0; i < D; i++)
            boundary[i] = aux::pair("-0.5","0.5");
    }

    void T(Vector &point) {
        point[2] += point[0];
        point[3] += point[1];

        Float coupling = xi*sin(2*aux::pi*(point[2] + point[3]));
        point[0] += k1*sin(2*aux::pi*point[2]) + coupling;
        point[1] += k2*sin(2*aux::pi*point[3]) + coupling;
        apply_boundary_conditions(point, bounding_box);
    }

    /*
    Mathematica code for computing the jacobian matrix:
     P1[p1_, p2_, q1_, q2_] := {
  p1 + K1/(2 \[Pi]) Sin[2 \[Pi] (q1 + p1)] +
   xi/(2 \[Pi]) Sin[2 \[Pi] (q1 + q2 + p1 + p2)],
  p2 + K2/(2 \[Pi]) Sin[2 \[Pi] (q2 + p2)] +
   xi/(2 \[Pi]) Sin[2 \[Pi] (q1 + q2 + p1 + p2)],
  q1 + p1,
  q2 + p2}

     D[P1[p1, p2, q1, q2], {{p1, p2, q1, q2}}] // MatrixForm
    */
    Matrix const& jacobian(Vector const& point) {
        _jacobian = Matrix::Zero(D,D);
        Float const& p1 = point[0];
        Float const& p2 = point[1];
        Float const& q1 = point[2];
        Float const& q2 = point[3];

        Float const& pi = aux::pi;
        Float coupling = 2*pi*xi*cos(2*pi*(p1 + q1 + p2 + q2));
        Float bla1 = 2*pi*k1*cos(2*pi*(p1 + q1));
        Float bla2 = 2*pi*k2*cos(2*pi*(p2 + q2));

        _jacobian(0,0) = 1 + bla1 + coupling;
        _jacobian(0,1) = coupling;
        _jacobian(0,2) = bla1 + coupling;
        _jacobian(0,3) = coupling;
        _jacobian(1,0) = coupling;
        _jacobian(1,1) = 1 + bla2 + coupling;
        _jacobian(1,2) = coupling;
        _jacobian(1,3) = bla2 + coupling;
        _jacobian(2,0) = _jacobian(2,2) = _jacobian(3,1) = _jacobian(3,3) = 1;
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (point[0] < "-0.4" or point[1] < "-0.4")
            return true;
        else
            return false;
    };
};


class Tent : public Map {
protected:
    Float a;
public:
    Tent(Float a) : Map(1, format("tent%.1f", a.toDouble())), a(a) {
        boundary[0] = aux::pair(0, 1);
    }

    void T(Vector & point) {
        if (point[0] < 1/a)
            point[0] *= a;
        else
            point[0] = a/(a - 1)*(1 - point[0]);
    }

    Matrix const& jacobian(Vector const& point) {
        if (point[0] < 1/a)
            _jacobian(0,0) = a;
        else
            _jacobian(0,0) = -a/(a - 1);
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (point[0] < "0.4")
            return true;
        else
            return false;
    };
};


//! Defined in Transient chaos: Complex dynamics in finite time scales (Lai + Tel)
class OpenTent : public Map {
    Float a;
    Float b;
public:
    OpenTent(Float a, Float b) : Map(1, format("tent%.f&%.f", a.toDouble(), b.toDouble())), a(a), b(b) {
        boundary[0] = aux::pair(0, 1);
    }

    void T(Vector & point) {
        if (point[0] < b/(a + b))
            point[0] *= a;
        else
            point[0] = b*(1 - point[0]);
    }

    Matrix const& jacobian(Vector const& point) {
        if (point[0] < b/(a + b))
            _jacobian(0,0) = a;
        else
            _jacobian(0,0) = b;
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (0 < point[0] and point[0] < 1)
            return false;
        else
            return true;
    };
};


class Logistic : public Map {
    Float r;
public:
    Logistic(Float r) : Map(1, format("logistic%1.f", r.toDouble())), r(r) {
        boundary[0] = aux::pair(0, 1);
    }

    void T(Vector & point) {
        point[0] = r*(1 - point[0])*point[0];
    }

    Matrix const& jacobian(Vector const& point) {
        _jacobian(0,0) = r*(1 - 2*point[0]);
        return _jacobian;
    }

    virtual bool has_exited(Vector const & point) {
        if (0 < point[0] and point[0] < Float("0.2"))
            return true;
        else
            return false;
    };
};


class NCoupledHenon : public Map {
protected:
    std::vector<Float> a;
    Float b;
    Float k;

public:
    NCoupledHenon(unsigned int D, Float min_a=3, Float max_a=5, Float b="0.3", Float k="0.4") : Map(D, format("ch%d", D)), k(k), b(b), a(D/2) {
        if(D % 2 != 0) {
            std::cout << "NCoupled dimension must be multiple of 2";
            exit(1);
        }
        // set values of a
        a[0] = min_a;
        a[D/2 - 1] = max_a;
        for (unsigned int i = 1; i < D/2 - 1; i++)
            a[i] = min_a + (max_a - min_a)*i/(D/2 - 1);

        for (unsigned int i = 0; i < D; i++)
            boundary[i] = aux::pair(-4, 4);
    }

    void T(Vector & point) {
        Float x0 = point[0];

        for (unsigned int i = 0; i < D/2; i++) {
            unsigned int iplus1 = (i + 1 + D/2)%(D/2);
            Float x = point[i];
            Float y = point[i + D/2];
            Float u = point[iplus1];
            if (iplus1 == 0)  // avoid retrieving already modified value
                u = x0;

            point[i] = a[i] - x*x + b*y;
            if (D > 2)
                point[i] += k*(x - u);
            point[i + D/2] = x;
        }
    }

    Matrix const& jacobian(Vector const& point) {
        _jacobian.setZero();

        for (unsigned int i = 0; i < D/2; i++) {
            unsigned int iplus1 = (i + 1 + D/2)%(D/2);

            _jacobian(i, i) = -2*point[i];
            _jacobian(i, i + D/2) = b;
            if (D > 2) {
                _jacobian(i, i) += k;
                _jacobian(i, iplus1) += -k;
            }
            _jacobian(i + D/2, i) = 1;
        }
        return _jacobian;
    }

    bool has_exited(Vector const& point) {
        for(unsigned int i = 0; i < D; i++)
            if (point[i] < -4 || point[i] > 4)
                return true;
        return false;
    }
};


}; // map

#endif /* defined(__chaospp__map__) */
