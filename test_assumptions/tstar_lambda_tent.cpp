/*

*/

#include <algorithm> // for max

#include "map.h"
#include "sampler.h"


double log_binomial_coefficient(int n, int k) {
    if (k == 0 or k == n)
        return 0;
    return n*log(n) - k*log(k) - (n - k)*log(n - k);
}

//! contains the entropy and lambda_L of the tent map
struct TentEntropy {
    unsigned int tobs;
    histogram::Histogram<double> histogram;
    std::vector<double> _entropy;
    double lambda_L;
    double const a;
    double const b;

    // "- 0001" avoids discretization errors because \lambdas are discrete.
    TentEntropy(unsigned int tobs, double a) : tobs(tobs), histogram(log(a/(a - 1)) - 0.0001, log(a) - 0.0001, tobs), _entropy(tobs + 1), a(a), b(a/(a - 1)) {
        for(unsigned int i = 0; i <= tobs; i++) {
            _entropy[i] = log_binomial_coefficient(tobs, i) + log(1/a)*i + log(1/b)*(tobs - i);
        }

        lambda_L = log(a)/a + log(b)/b;
    }

    double entropy(double value) {
        return _entropy[histogram.bin(value)];
    }

    double d_entropy(double value) const {
        double i = histogram.bin(value);

        if (i == 0 or i == tobs)
            return 0;
        // the analytical derivative of the equation in `_entropy`
        return log(1/a) - log(1/b) - log(i) + log(tobs - i);
    }
};


//! the proposal distribution given by the tstar from the thesis
template <typename Observable>
class TentTstarProposal : public proposal::Isotropic<Observable> {
    Float delta_0;
    unsigned int tobs;

public:
    TentEntropy ent;

    // according to formula in the thesis.
    Float sigma(Observable const& result) const {
        double lambda = result.lyapunov();

        // tent properties
        double lambda_L = this->ent.lambda_L;
        double S_derivative = this->ent.d_entropy(lambda);

        //S_derivative = 1; // beta = 1
        double delta_t = 1/fabs(S_derivative*(lambda - lambda_L));

        double t_star = std::max(0.0, tobs - delta_t);

        return delta_0*exp(-lambda*t_star);
    }

    TentTstarProposal(unsigned int tobs, std::vector<aux::pair> const& boundary, Float delta_0, double a) : tobs(tobs),
            proposal::Isotropic<Observable>(boundary), delta_0(delta_0), ent(tobs, a) {}
};


void measure(unsigned int tobs, double delta, std::string directory="./") {
    mpfr::mpreal::set_default_prec(64);

    std::string map_name("tent3");

    map::Tent map(3);

    TentTstarProposal<observable::Lyapunov> proposal(tobs, map.boundary, delta, 3);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < 100000; i++) {
        Vector state = proposal.proposeUniform();

        observable::Lyapunov result(map, tobs);
        result.observe(state);

        Vector state_prime = proposal.propose(result);

        observable::Lyapunov result_prime(map, tobs);
        result_prime.observe(state_prime);

        std::vector<double> temp = {result.state[0].toDouble(), result.lyapunov(), proposal.ent.entropy(result.lyapunov()),
                                    result_prime.state[0].toDouble(), result_prime.lyapunov(), proposal.ent.entropy(result_prime.lyapunov()),
                                    proposal.log_acceptance(result, result_prime)};
        data.push_back(temp);
    }

    io::save(data, directory + map_name + format("_%d_fh_%.2f", tobs, delta) + ".dat");
}


int main() {
    measure(10, 0.1);
    measure(20, 0.1);
    measure(30, 0.1);

    measure(10, 1);
    measure(20, 1);
    measure(30, 1);

    measure(10, 10);
    measure(20, 10);
    measure(30, 10);
}
