#ifndef chaospp_tstar_proposal_h
#define chaospp_tstar_proposal_h

#include <vector>

double finite_difference(std::vector<double> const& values, unsigned int index) {
    if (index == values.size() - 1)
        return values[index] - values[index - 1];
    else if (index == 0)
        return values[index + 1] - values[index];
    else
        return (values[index + 1] - values[index - 1])/2;
}


//! the proposal distribution given by the tstar from the thesis
template <typename Observable>
class TstarProposal : public proposal::Isotropic<Observable> {
    SamplingHistogram<Observable> const& histogram;
    Float delta_0;
    unsigned int tobs;

    double t_star(Observable const& result) const {
        double lambda = result.lyapunov();
        unsigned int bin = histogram.bin(result.observable());

        // lambda that maximizes \pi(E)
        unsigned int bin_max = 0;
        double max = -1000000;
        for (unsigned int b = 0; b <= histogram.bins(); b++) {
            if (histogram.entropy(b) > max) {
                bin_max = b;
                max = histogram.entropy(b);
            }
        }

        double lambda_L = histogram.value(bin_max);

        // compute derivative
        double d_log_pi;
        if (bin == histogram.log_pi.size() - 1)
            d_log_pi = histogram.log_pi[bin] - histogram.log_pi[bin - 1];
        else if (bin == 0)
            d_log_pi = histogram.log_pi[bin + 1] - histogram.log_pi[bin];
        else
            d_log_pi = (histogram.log_pi[bin + 1] - histogram.log_pi[bin - 1])/2;
        d_log_pi /= this->histogram.h();

        double delta_t = 1/fabs(d_log_pi*(lambda - lambda_L));

        if (std::isnan(d_log_pi) or std::isinf(d_log_pi))
            delta_t = tobs;

        return std::max(0.0, tobs - delta_t);
    }

    // according to formula in the thesis.
    Float sigma(Observable const& result) const {
        double lambda = result.lyapunov();

        return delta_0*exp(-lambda*this->t_star(result));
    }

public:
    TstarProposal(std::vector<aux::pair> const& boundary, Float delta_0, unsigned int tobs, SamplingHistogram<Observable> const& histogram) :
            proposal::Isotropic<Observable>(boundary), delta_0(delta_0), tobs(tobs), histogram(histogram) {}
};


#endif
