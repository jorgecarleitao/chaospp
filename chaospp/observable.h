#ifndef chaospp_observable_h
#define chaospp_observable_h

#include "auxiliar.h"
#include "map.h"
#include <Eigen/Eigenvalues>


class ComputeMatrix {
public:
    Matrix jacobian;  //! The final jacobian matrix.

    ComputeMatrix(map::Map & map) : jacobian(Matrix::Identity(map.D, map.D)) {}

    void evolve(map::Map & map, Vector & point) {
        jacobian *= map.jacobian(point);
    }

    Float stretch() const {
        return abs(es.eigenvalues()[i_max]);
    }

    Vector eigenvector() const {
        Vector vector(jacobian.rows());
        for (unsigned int i = 0; i < jacobian.rows(); i++)
            vector[i] = es.eigenvectors().col(i_max)[i];
        return vector;
    }

protected:
    unsigned int i_max;
    Eigen::EigenSolver<Matrix> es;

    void finalise(map::Map & map) {
        es.compute(jacobian);

        Float maximum(-1000);
        for (unsigned int i = 0; i < map.D; i++)
            if (abs(es.eigenvalues()[i]) > maximum) {
                i_max = i;
                maximum = abs(es.eigenvalues()[i]);
            }
    }

    void initialize() {
        jacobian = Matrix::Identity(jacobian.rows(), jacobian.rows());
    }

    ComputeMatrix & operator=(ComputeMatrix const& other) {
        this->jacobian = other.jacobian;
        this->i_max = other.i_max;
        this->es = other.es;
        return *this;
    }
};


namespace observable {

//! The outcome of the evolution of the system. This class calls map iterations and stores relevant intermediate results.
//! It contains a single attribute, `state`, the initial state.
template <typename T>
class Observable {
public:
    Vector state; // the initial state. Set in "observe".
    typedef T Type;

    //! The initial state of the map that fully characterizes the system.
    virtual void observe(Vector const& state) {
        this->state = state;
    }

    virtual Observable & operator=(Observable const& other) {
        this->state = other.state;
        return *this;
    }

    virtual T observable() const = 0;
};


//! the escape time of the state for open systems.
class EscapeTime : public Observable<unsigned int> {
protected:
    virtual void finalize() {}

    virtual void initialize() {
        escape_time = 0;
    }

public:
    map::Map & map;

    //! The escape time of `state`. It is computed on "observe"
    unsigned int escape_time;
    unsigned int max_time;

    EscapeTime(map::Map & map, unsigned int max_time=std::numeric_limits<unsigned int>::max()) : map(map), max_time(max_time), escape_time(0) {}

    virtual bool has_exited(Vector const& point) const {return map.has_exited(point);}

    //! Escape time function
    void observe(Vector const& state) {
        Observable::observe(state);
        initialize();

        Vector point = state;
        evolve(point);
        while (not has_exited(point) and escape_time < max_time) {
            evolve(point);
        }
        finalize();
    }

    virtual void evolve(Vector & point) {
        map.T(point);
        escape_time++;
    }

    EscapeTime & operator=(EscapeTime const& other) {
        Observable::operator=(other);
        this->escape_time = other.escape_time;
        this->max_time = other.max_time;
        return *this;
    }

    virtual unsigned int observable() const {
        return escape_time;
    }
};

const Vector UNDEFINED_VECTOR;

//! Computes the escape time of the state and its FT Lyapunov exponent
//! It is constructed from a map, an optional max_time (default: infinite), and an optional tangent vector (default: random vector).
//! It evolves the system and tangent vector in time until the state escapes or up to max_time.
class EscapeWithVector : public EscapeTime {
protected:
    //! The tangent vector after `escape_time` steps.
    Vector tangent;

    virtual void initialize() {
        EscapeTime::initialize();
        tangent = aux::unitaryVector(map.D);
    }

public:

    EscapeWithVector(map::Map & map, unsigned int max_time=std::numeric_limits<unsigned int>::max(), Vector const& tangent=UNDEFINED_VECTOR) :
            EscapeTime(map, max_time), tangent(tangent) {
        if (this->tangent == UNDEFINED_VECTOR)
            this->tangent = aux::unitaryVector(map.D);
    }

    virtual Float stretch() const {
        return aux::get_norm(tangent);
    }

    double lyapunov() const {
        return ((double) log(stretch()))/this->escape_time;
    }

    virtual void evolve(Vector & point) {
        map.dT(point, tangent);
        EscapeTime::evolve(point);
    }

    EscapeWithVector & operator=(EscapeWithVector const& other) {
        EscapeTime::operator=(other);
        this->tangent = other.tangent;
        return *this;
    }
};


//! Computes the escape time of the state and the Jacobian matrix of the trajectory.
class EscapeWithMatrix : public EscapeTime, public ComputeMatrix {
public:
    EscapeWithMatrix(map::Map & map, unsigned int max_time=std::numeric_limits<unsigned int>::max()) : EscapeTime(map, max_time), ComputeMatrix(map) {}

    double lyapunov() const {
        return log(stretch()).toDouble()/this->escape_time;
    }

    void evolve(Vector & point) {
        ComputeMatrix::evolve(map, point);
        EscapeTime::evolve(point);
    }

    void finalize() {
        ComputeMatrix::finalise(map);
        EscapeTime::finalize();
    }

    EscapeWithMatrix & operator=(EscapeWithMatrix const& other) {
        EscapeTime::operator=(other);
        ComputeMatrix::operator=(other);
        return *this;
    }
};


//! Computes the finite time Lyapunov exponent
class Lyapunov : public Observable<double>, public ComputeMatrix {
public:
    map::Map & map;
    unsigned int tobs;

    Lyapunov(map::Map & map, unsigned int tobs) : map(map), ComputeMatrix(map), tobs(tobs) {}

    virtual void finalise(Vector const&) {
        ComputeMatrix::finalise(map);
    }

    //! Lyapunov exponent function
    virtual void observe(Vector const& state) {
        Observable::observe(state);
        ComputeMatrix::initialize();

        Vector point = state;

        for(unsigned int escape_time = 0; escape_time < tobs; escape_time++) {
            ComputeMatrix::evolve(map, point);
            map.T(point);
        }

        finalise(point);
    }

    double lyapunov() const {
        return log(stretch()).toDouble()/tobs;
    }

    Lyapunov & operator=(Lyapunov const& other) {
        Observable::operator=(other);
        ComputeMatrix::operator=(other);
        tobs = other.tobs;
        return *this;
    }

    virtual double observable() const {
        return lyapunov();
    }
};

}

#endif
