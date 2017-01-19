#ifndef chaospp_test_isotropic_proposal_h
#define chaospp_test_isotropic_proposal_h

#include "proposal.h"


class State : public observable::Observable<Vector> {
    virtual Vector observable() const {
        return state;
    }
};


class ConstantSigmaProposal : public proposal::Isotropic<State> {
public:
    Float sigma0;
    ConstantSigmaProposal(Float sigma0, std::vector<aux::pair> const& boundary) :
            proposal::Isotropic<State>(boundary), sigma0(sigma0) {}

    Float sigma(State const&) const {
        return sigma0;
    }
};


TEST(IsotropicProposal, 1D) {
    // Test that
    std::vector<aux::pair> boundary(1);
    boundary[0] = aux::pair(-10, 10);

    ConstantSigmaProposal proposal(0.2, boundary);

    Vector point(1), point_prime;
    point[0] = 0;

    State r;
    r.observe(point);

    double mean = 0;
    const unsigned int samples = 1000;
    for(unsigned int i = 0; i < samples; i++) {
        point_prime = proposal.propose(r);
        // compute mean of |x - x'|
        double value = fabs(proposal.get_delta().toDouble());
        mean += value;
    }
    mean /= samples;

    // expected: mean = 0.2 within 3 sigmas
    EXPECT_NEAR(proposal.sigma0.toDouble(), mean, 3*proposal.sigma0.toDouble()/sqrt(samples));
}


TEST(IsotropicProposal, 2D) {
    // Test that
    static const unsigned int D = 2;
    std::vector<aux::pair> boundary(2);
    boundary[0] = aux::pair(-10, 10);
    boundary[1] = aux::pair(-10, 10);

    ConstantSigmaProposal proposal(0.2, boundary);

    Vector point(2), point_prime(2);
    point[0] = 0;
    point[1] = 0;

    State r;
    r.observe(point);

    double mean = 0;
    const unsigned int samples = 1000;
    for(unsigned int i = 0; i < samples; i++) {
        point_prime = proposal.propose(r);
        // compute mean of |x - x'|
        double value = fabs(proposal.get_delta().toDouble());
        mean += value;
    }
    mean /= samples;

    // expected: mean = 0.2 within 3 sigmas
    EXPECT_NEAR(proposal.sigma0.toDouble(), mean, 3*proposal.sigma0.toDouble()/sqrt(samples));
}


class SigmaProposal : public proposal::Isotropic<State> {
public:
    SigmaProposal(std::vector<aux::pair> const& boundary) : proposal::Isotropic<State>(boundary) {}

    Float sigma(State const& result) const {
        if (result.state[0] == 0)
            return 10;
        else
            return 5;
    }
};


TEST(Acceptance, 2D) {
    std::vector<aux::pair> boundary(2);
    boundary[0] = aux::pair(-10, 10);
    boundary[1] = aux::pair(-10, 10);

    SigmaProposal proposal(boundary);

    Vector point(2), point_prime(2);
    point[0] = 0;
    point[1] = 0;

    State r, rP;
    r.observe(point);

    double mean = 0, meanD = 0, meanD2 = 0;
    const unsigned int samples = 10;
    for(unsigned int i = 0; i < samples; i++) {
        point_prime = proposal.propose(r);
        rP.observe(point_prime);
        // compute mean of |x - x'|
        double delta = proposal.get_delta().toDouble();

        mean += proposal.log_acceptance(r, rP);
        meanD += fabs(delta);
        meanD2 += delta*delta/(aux::pi/2).toDouble();
    }
    mean /= samples;
    meanD /= samples;
    meanD2 /= samples;

    EXPECT_NEAR(log(2) - 0.5*(4 - 1)*meanD2/100, mean, 0.1);
    EXPECT_NEAR(10, meanD, 3*10/sqrt(samples));
}


#endif
