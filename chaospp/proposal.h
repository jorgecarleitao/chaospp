#ifndef chaospp_proposal_h
#define chaospp_proposal_h

#include "assert.h"

#include "auxiliar.h"
#include "observable.h"

namespace proposal {

//! sends the point to inside the boundary
void bound_initial_condition(Vector & point, std::vector<aux::pair> const& boundary) {
    assert(boundary.size() == point.size());
    for(unsigned int i = 0; i < point.size(); i++) {
        while(point[i] > boundary[i].second)
            point[i] -= boundary[i].second - boundary[i].first;
        while(point[i] < boundary[i].first)
            point[i] += boundary[i].second - boundary[i].first;
    }
}

Vector proposeUniform(std::vector<aux::pair> const& boundary) {
    Vector proposal(boundary.size());
    for(unsigned int i = 0; i < boundary.size(); i++) {
        aux::pair const& box = boundary[i];
        proposal[i] = box.first + (box.second - box.first)*aux::urandom();
    }
    return proposal;
}

Vector proposeIsotropic(Vector point, Vector const& vector, Float const& sigma, std::vector<aux::pair> const& boundary) {
    for(unsigned d = 0; d < point.size(); d++)
        point[d] += sigma*vector[d];

    bound_initial_condition(point, boundary);
    return point;
}

inline double logAcceptanceIsotropic(Float const& sigma, Float const& sigmaPrime, Float delta) {
    double ratio = Float(delta/sigma).toDouble();
    double ratio_sigma = Float(sigma/sigmaPrime).toDouble();

    return log(ratio_sigma) - 0.5*ratio*ratio*(ratio_sigma*ratio_sigma - 1);
}

Vector proposeAnisotropic(Vector point, Matrix const& jacobian, Float const& sigma0, std::vector<aux::pair> const& boundary) {
    Eigen::JacobiSVD<Matrix> svd(jacobian, Eigen::ComputeFullV);

    Matrix v_matrix = svd.matrixV();
    Vector const& singular_values = svd.singularValues();

    Vector delta = aux::unitaryVector((unsigned int)point.size());

    // build vector
    // build (Sigma^+)^-1*vector
    for (unsigned int d = 0; d < point.size(); d++) {
        if(singular_values(d) > 1)
            delta(d) *= sigma0/singular_values(d);
        else
            delta(d) = 0;
    }

    // V*(Sigma^+)^-1*vector
    delta = v_matrix*delta;

    point += delta;

    bound_initial_condition(point, boundary);
    return point;
}


// A generic class that implements proposals.
template <typename Observable>
class Proposal {
    static_assert(
            std::is_base_of<observable::Observable<typename Observable::Type>, Observable>::value,
            "Observable must be a subclass of observable::Observable"
    );
protected:
    Float delta;
    std::vector<aux::pair> const& boundary;
    unsigned int D;
public:

    Proposal(std::vector<aux::pair> const& boundary) : boundary(boundary), D((unsigned int) boundary.size()) {}

    Vector proposeUniform() {
        return proposal::proposeUniform(boundary);
    }

    virtual Vector propose(Observable const& result) = 0;

    virtual double log_acceptance(Observable const&, Observable const&) const = 0;

    virtual void update(Observable const&, Observable const&) {}

    Float const& get_delta() const {
        return delta;
    }
};


// Uniform proposal on the boundary region
template <typename Observable>
class Uniform : public Proposal<Observable> {
public:
    Uniform(std::vector<aux::pair> const& boundary) : Proposal<Observable>(boundary) {}

    Vector propose(Observable const& result) {
        Vector newState = this->proposeUniform();
        this->delta = aux::get_norm(newState - result.state);
        return newState;
    }

    double log_acceptance(Observable const&, Observable const&) const {
        return 0;
    }
};


// "Power Law isotropic proposal"
// aka "Exponential Stagger Distribution" in http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.2261
template <typename Observable>
class PowerLawIsotropic : public Proposal<Observable> {
protected:
    Float delta;
    Float min_s, max_s;
public:

    PowerLawIsotropic(std::vector<aux::pair> const& boundary, Float min_s, Float max_s) :
            Proposal<Observable>(boundary), min_s(-min_s), max_s(-max_s) {}

    virtual Vector propose(Observable const& result) {
        delta = exp(min_s + (max_s - min_s)*aux::urandom());
        return proposeIsotropic(result.state, aux::unitaryVector(this->D), delta, this->boundary);
    }

    virtual double log_acceptance(Observable const&, Observable const&) const {
        return 0;
    }
};


// Half-normal isotropic proposal. Overload "sigma(observable)" to define the function
template <typename Observable>
class Isotropic : public Proposal<Observable> {
public:
    Isotropic(std::vector<aux::pair> const& boundary) : Proposal<Observable>(boundary) {}

    virtual Float sigma(Observable const& result) const = 0;

    virtual Vector propose(Observable const& result) {
        static const Float constant = sqrt(aux::pi/2);
        // we multiply here by constant, and divide in the acceptance respectively
        this->delta = sigma(result)*constant*abs(aux::nrandom());
        return proposeIsotropic(result.state, aux::unitaryVector(this->D), this->delta, this->boundary);
    }

    virtual double log_acceptance(Observable const& result, Observable const& result_prime) const {
        // we divide here by constant, and multiply in propose respectively
        static const Float constant = sqrt(aux::pi/2);
        return logAcceptanceIsotropic(sigma(result), sigma(result_prime), this->delta/constant);
    }
};


// Lyapunov proposal of http://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.052916
template <typename Observable>
class LyapunovIsotropic : public Isotropic<Observable> {
protected:
    Float sigma0;
public:

    virtual Float sigma(Observable const& result) const {
        return sigma0/result.stretch();
    }

    LyapunovIsotropic(std::vector<aux::pair> const& boundary, Float const& sigma0=10) :
            Isotropic<Observable>(boundary), sigma0(sigma0) {}
};


//! "Adaptive proposal" of http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.220601
template <typename Observable=observable::EscapeTime>
class Adaptive : public Isotropic<Observable> {
    static_assert(
            std::is_base_of<observable::EscapeTime, Observable>::value,
            "Observable must be a subclass of EscapeTime"
    );

protected:
    Float _sigma;
    Float factor;
    Float max_sigma;
public:
    Adaptive(std::vector<aux::pair> const& boundary, Float const& factor = "1.1") :
        Isotropic<Observable>(boundary), _sigma(1), factor(factor), max_sigma(10) {}

    virtual Float sigma(Observable const& result) const {
        return _sigma;
    }

    virtual void update(Observable const& result, Observable const& result_prime) {
        if (result_prime.escape_time >= result.escape_time) {
            if (_sigma * factor > max_sigma)
                _sigma = max_sigma;
            else
                _sigma *= factor;
        }
        else
            _sigma /= factor;
    }
};


// Anisotropic proposal such that the proposal is isotropic in the end of the trajectory
template <typename Observable=observable::EscapeWithMatrix>
class Anisotropic : public Proposal<Observable> {
public:
    Anisotropic(std::vector<aux::pair> const& boundary) : Proposal<Observable>(boundary) {}

    virtual Vector propose(Observable const& result) {
        return proposeAnisotropic(result.state, result.jacobian, 10, this->boundary);
    };

    // We have not done the calculation for this acceptance.
    // todo: add formula here.
    virtual double log_acceptance(Observable const&, Observable const&) const {
        assert(1 == 0);
        return 0;
    }
};

}

#endif
