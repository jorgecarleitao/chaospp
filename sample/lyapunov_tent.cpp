/*
 * Measures the finite-time Lyapunov exponent.
 */
#include <algorithm> // for max

#include "map.h"
#include "sampler.h"
#include "observable.h"
#include "thesis_proposal.h"


double log_binomial_coefficient(int n, int k) {
    if (k == 0 or k == n)
        return 0;
    return n*log(n) - k*log(k) - (n - k)*log(n - k);
}


class TestHistogram : public SamplingHistogram<observable::Lyapunov> {
    typedef observable::Lyapunov Observable;
protected:
    std::vector<double> acceptance;
    std::vector<double> acceptance2;

    std::vector<double> delta_e;
    std::vector<double> delta_e2;
    std::vector<double> conditional_histogram;
    std::vector<double> log_delta_histogram;
    std::vector<double> log_delta2_histogram;
    std::vector<double> lambda_histogram;
    std::vector<double> lambda2_histogram;

    void export_acceptance(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (count > 0) {
                std::vector<double> row(3);
                row[0] = this->value(bin);
                row[1] = acceptance[bin]*1./count;

                row[2] = acceptance2[bin]*1./count;
                row[2] = row[2] - row[1]*row[1];
                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

    void export_delta_e(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (count > 0) {
                std::vector<double> row(3);
                row[0] = this->value(bin);
                row[1] = delta_e[bin]*1./count;

                row[2] = delta_e2[bin]*1./count;
                row[2] = row[2] - row[1]*row[1];
                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

    void export_conditional(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (count > 0) {
                std::vector<double> row(2);
                row[0] = this->value(bin);
                row[1] = conditional_histogram[bin]*1./count;
                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

    void export_log_delta(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (this->conditional_histogram[bin] > 0) {
                std::vector<double> row(5);
                row[0] = this->value(bin);
                row[1] = log_delta_histogram[bin]*1./this->conditional_histogram[bin];
                // variance
                row[2] = log_delta2_histogram[bin]*1./this->conditional_histogram[bin];
                row[2] -= row[1]*row[1];
                // std
                row[2] = sqrt(row[2]);

                row[3] = lambda_histogram[bin]*1./count;
                row[4] = lambda2_histogram[bin]*1./count;
                // std
                row[4] -= row[3]*row[3];
                row[4] = sqrt(row[4]);
                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

public:

    TestHistogram(double lowerBound, double upperBound, unsigned int bins)
            : SamplingHistogram<observable::Lyapunov>(lowerBound, upperBound, bins),
              acceptance(this->bins()),
              acceptance2(this->bins()),
              delta_e(this->bins()),
              delta_e2(this->bins()),
              log_delta_histogram(this->bins()),
              log_delta2_histogram(this->bins()),
              conditional_histogram(this->bins()),
              lambda_histogram(this->bins()),
              lambda2_histogram(this->bins())
    {}

    virtual void measure(Observable const& result, Observable const& result_prime, double accept) {
        SamplingHistogram<observable::Lyapunov>::measure(result, result_prime, accept);
        unsigned int bin = this->bin(result.observable());

        acceptance[bin] += accept;
        acceptance2[bin] += accept*accept;

        double diff = result.tobs*(result_prime.lyapunov() - result.lyapunov());
        delta_e[bin] += diff;
        delta_e2[bin] += diff*diff;

        if (0.1 < accept and accept < 0.9) {
            conditional_histogram[bin] += 1;
            Float norm = aux::get_norm(result_prime.state - result.state);
            double value = -log(norm).toDouble();
            log_delta_histogram[bin] += value;
            log_delta2_histogram[bin] += value * value;
        }
        double value = log(result.stretch()).toDouble();
        lambda_histogram[bin] += value;
        lambda2_histogram[bin] += value * value;
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        SamplingHistogram<observable::Lyapunov>::export_pretty(file_name, directory);
        export_acceptance("acceptance_" + file_name, directory);
        export_delta_e("deltae_" + file_name, directory);
        export_log_delta("log_delta_" + file_name, directory);
        export_conditional("conditional_acceptance_" + file_name, directory);
    }

    void reset() {
        SamplingHistogram<observable::Lyapunov>::reset();

        // reset to 0
        std::fill(acceptance.begin(), acceptance.end(), 0);
        std::fill(acceptance2.begin(), acceptance2.end(), 0);
        std::fill(delta_e.begin(), delta_e.end(), 0);
        std::fill(delta_e2.begin(), delta_e2.end(), 0);
        std::fill(log_delta_histogram.begin(), log_delta_histogram.end(), 0);
        std::fill(log_delta2_histogram.begin(), log_delta2_histogram.end(), 0);
        std::fill(conditional_histogram.begin(), conditional_histogram.end(), 0);
        std::fill(lambda_histogram.begin(), lambda_histogram.end(), 0);
        std::fill(lambda2_histogram.begin(), lambda2_histogram.end(), 0);
    }
};


// Thesis Proposal
void fh_tp(unsigned int tobs, double delta, std::string directory="./") {

    mpfr::mpreal::set_default_prec(128);

    double a = 3;
    map::Tent map(a);
    observable::Lyapunov observable(map, tobs);

    TestHistogram histogram(log(a / (a - 1)) - 0.0001, log(a) - 0.0001, tobs);

    TstarProposal<observable::Lyapunov> proposal(map.boundary, delta, tobs, histogram);

    std::vector<double> entropy(tobs + 1);
    for(unsigned int i = 0; i <= tobs; i++) {
        entropy[i] = log_binomial_coefficient(tobs, i) + log(1/a)*i + log(1/(a / (a - 1)))*(tobs - i);
    }
    histogram.set_entropy(entropy);

    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -histogram.entropy(bin);
    }

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(100000);

    histogram.export_histogram(format("tent3_%d_fh_tp1_%.2f.dat", tobs, delta), directory);
}


// Uniform Proposal
void fh_up(unsigned int tobs, std::string directory="./") {

    mpfr::mpreal::set_default_prec(128);

    double a = 3;
    map::Tent map(a);
    observable::Lyapunov observable(map, tobs);

    TestHistogram histogram(log(a / (a - 1)) - 0.0001, log(a) - 0.0001, tobs);

    proposal::Uniform<observable::Lyapunov> proposal(map.boundary);

    std::vector<double> entropy(tobs + 1);
    for(unsigned int i = 0; i <= tobs; i++) {
        entropy[i] = log_binomial_coefficient(tobs, i) + log(1/a)*i + log(1/(a / (a - 1)))*(tobs - i);
    }
    histogram.set_entropy(entropy);

    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -histogram.entropy(bin);
    }

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(10000);

    histogram.export_histogram(format("tent3_%d_fh_up.dat", tobs), directory);
}


// Thesis Proposal in uniform sampling
void us_tp(unsigned int tobs, double delta, std::string directory="./") {

    mpfr::mpreal::set_default_prec(128);

    double a = 3;
    map::Tent map(a);
    observable::Lyapunov observable(map, tobs);

    TestHistogram histogram(log(a / (a - 1)) - 0.0001, log(a) - 0.0001, tobs);

    std::vector<double> entropy(tobs + 1);
    for(unsigned int i = 0; i <= tobs; i++) {
        entropy[i] = log_binomial_coefficient(tobs, i) + log(1/a)*i + log(1/(a / (a - 1)))*(tobs - i);
    }
    histogram.set_entropy(entropy);

    TstarProposal<observable::Lyapunov> proposal(map.boundary, delta, tobs, histogram);

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(10000);

    histogram.export_histogram(format("tent3_%d_us_tp_%.2f.dat", tobs, delta), directory);
}


//// Power-law proposal
void fh_sy(unsigned int tobs, std::string directory="./") {

    mpfr::mpreal::set_default_prec(128);

    double a = 3;
    map::Tent map(a);

    observable::Lyapunov observable(map, tobs);

    proposal::PowerLawIsotropic<observable::Lyapunov> proposal(map.boundary, -3, 40);

    TestHistogram histogram(log(a/(a - 1)) - 0.0001, log(a) - 0.0001, tobs);

    std::vector<double> entropy(tobs + 1);
    for(unsigned int i = 0; i <= tobs; i++) {
        entropy[i] = log_binomial_coefficient(tobs, i) + log(1/a)*i + log(1/(a / (a - 1)))*(tobs - i);
    }
    histogram.set_entropy(entropy);

    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -histogram.entropy(bin);
    }

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(100000);

    histogram.export_histogram(format("tent3_%d_fh_sy.dat", tobs), directory);
}


int main() {
//    fh_tp(10, 1);
//    fh_tp(20, 1);
//    fh_tp(30, 1);

//    fh_up(10);
//    fh_up(20);
//    fh_up(30);

//    us_tp(10, 1);
//    us_tp(20, 1);
//    us_tp(30, 1);

//    us_tp(20, 10);

    //fh_sy(10);
    //fh_sy(20);
    fh_sy(30);
    return 0;
}
