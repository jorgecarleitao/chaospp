#ifndef chaospp_sampler_h
#define chaospp_sampler_h

#include <algorithm>

#include "proposal.h"
#include "map.h"
#include "histogram.h"


//! This is an histogram that contains
template <typename Observable>
class SamplingHistogram : public histogram::Histogram<typename Observable::Type> {
    typedef typename Observable::Type T;
protected:
    bool has_exact_entropy;  // the main histogram
    std::vector<double> _entropy; // in case has_exact_entropy
public:
    std::vector<double> log_pi;  // log of the sampling distribution

    SamplingHistogram(T lowerBound, T upperBound, unsigned int bins) :
            histogram::Histogram<T>(lowerBound, upperBound, bins), log_pi(bins + 1, 0), _entropy(bins + 1),
            has_exact_entropy(false) {}

    virtual void measure(Observable const& result, Observable const&, double) {
        this->add(result.observable());
    }

    using histogram::Histogram<T>::bin;

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        histogram::Histogram<T>::export_histogram("histogram_" + file_name, directory);
    }

    virtual void export_pretty(std::string file_name, std::string directory="") const {
        histogram::Histogram<T>::export_pretty("histogram_" + file_name, directory);
    }

    double entropy(unsigned int b) const {
        if (has_exact_entropy)
            return _entropy[b];
        else
            return log((*this)[b]) - log_pi[b];
    }

    void set_entropy(std::vector<double> const& entropy) {
        assert(entropy.size() == this->_entropy.size());
        this->_entropy = entropy;
        has_exact_entropy = true;
    }

    //! exports the best estimator of the normalized entropy, S(E) : \sum(\exp(S(E))) == 1
    void export_entropy(std::string file_name, std::string directory="") const {
        std::vector<std::vector<double> > data;

        // P(E) = H(E)/pi(E)*exp(-C) = 1  (C is the constant)
        // \sum[exp(log_H - log_pi)] = exp(C)
        // \sum[exp(log_H - log_pi + a_max - a_max)] = exp(C)
        // C = a_max + log(\sum[\exp(log_H - log_pi - a_max])

        // compute the maximum of log_H - log_pi, a_max:
        double a_max = -100;
        for (unsigned int b = 0; b <= this->bins(); b++) {
            double a = entropy(b);
            if (a > a_max)
                a_max = a;
        }

        // compute the normalization C:
        double sum = 0;
        for (unsigned int b = 0; b <= this->bins(); b++) {
            sum += exp(entropy(b) - a_max);
        }
        double C = a_max + log(sum);

        // export the entropy:
        for (unsigned int b = 0; b <= this->bins(); b++) {
            std::vector<double> row(2);
            row[0] = this->value(b);
            row[1] = entropy(b) - C;
            data.push_back(row);
        }
        io::save(data, directory + "entropy_" + file_name);
    }
};


//! This is a general class that implements Metropolis-Hastings.
//! It requires the observable over which the algorithm is going to be used,
//! a proposal distribution,
//! and an histogram that discretizes the observable and defines the sampling distribution.
template <typename Observable>
class MetropolisHastings {
protected:
    typedef SamplingHistogram<Observable> Histogram;
    typedef proposal::Proposal<Observable> Proposal;

    Observable const& observable;
    Proposal & proposal;
    Histogram & histogram;

    //! returns log(pi'/pi) + log(g'/g)
    double log_acceptance(Observable const& result, Observable const& result_prime) const {
        unsigned int bin = histogram.bin(result.observable());
        unsigned int bin_prime = histogram.bin(result_prime.observable());

        double delta = histogram.log_pi[bin_prime] - histogram.log_pi[bin];
        return delta + this->proposal.log_acceptance(result, result_prime);
    }

public:

    MetropolisHastings(Observable const& observable, Proposal & proposal, Histogram & histogram) :
            observable(observable), proposal(proposal), histogram(histogram) {}

    virtual void measure(Observable const& result, Observable const& result_prime, double acceptance) {
        histogram.measure(result, result_prime, acceptance);
    }

    inline Observable propose(Observable & result) {
        // generate point x' and observables E'
        Observable result_prime(result);
        result_prime.observe(proposal.propose(result));

        proposal.update(result, result_prime);

        while (histogram.invalid_value(result_prime.observable())) {
            result_prime.observe(proposal.propose(result));
        }
        return result_prime;
    }

    void markov_step(Observable & result, bool measure=true) {
        // generate proposal
        Observable result_prime(this->propose(result));

        // compute acceptance
        double log_acceptance = this->log_acceptance(result, result_prime);
        double acceptance = std::min(1.0, exp(log_acceptance));

        // measure previous state (result), proposed state (result_prime) and acceptance
        if (measure)
            this->measure(result, result_prime, acceptance);

        // accept/reject
        if (aux::urandom() < acceptance)
            result = result_prime;
    }

    //! performs a round-trip, from minBin to maxBin.
    void round_trip(Observable & result, unsigned int minBin=1, unsigned int maxBin=0) {
        bool goingUp = false;

        if (maxBin == 0)
            maxBin = histogram.bins() - 1;

        while (true) {
            this->markov_step(result);

            unsigned int bin = histogram.bin(result.observable());

            if (!goingUp and bin == maxBin)
                goingUp = true;
            if (goingUp and bin == minBin)
                break;
        }
    }

    virtual void sample(unsigned int total_samples, unsigned int convergence_samples=0) {
        Observable result(observable);
        result.observe(this->proposal.proposeUniform());

        // reach assymptotic distribution
        for(unsigned int sample = 0; sample < convergence_samples; sample++)
            this->markov_step(result, false);

        // sample
        for (unsigned int sample = 0; sample < total_samples; sample++) {
            this->markov_step(result);
        }
    }
};


template <typename Observable>
class WangLandau : public MetropolisHastings<Observable> {
protected:
    typedef SamplingHistogram<Observable> Histogram;
    typedef proposal::Proposal<Observable> Proposal;

    double f;
public:

    WangLandau(Observable const& observable, Proposal & proposal, Histogram & histogram) :
            MetropolisHastings<Observable>(observable, proposal, histogram), f(1) {}

    virtual void measure(Observable const& result, Observable const& result_prime, double acceptance) {
        MetropolisHastings<Observable>::measure(result, result_prime, acceptance);

        unsigned int bin = this->histogram.bin(result.observable());
        this->histogram.log_pi[bin] -= f; // Wang-Landau step (S+=f <=> log_pi-=f)
    }

    void sample(unsigned int steps, unsigned int total_samples) {
        Observable result(this->observable);
        result.observe(this->proposal.proposeUniform());

        for (unsigned int step = 0; step < steps; step++) {
            this->histogram.reset();
            for (unsigned int sample = 0; sample < total_samples; sample++) {
                this->markov_step(result);
            }
            f /= 2;
        }
    }

    void approximate_entropy(unsigned int steps, unsigned int round_trips) {
        Observable result(this->observable);
        result.observe(this->proposal.proposeUniform());

        for (unsigned int step = 0; step < steps; step++) {
            this->histogram.reset();
            for (unsigned int round_trip = 0; round_trip < round_trips; round_trip++) {
                std::cout << format("%d/%d", round_trip, round_trips) << std::endl;
                this->round_trip(result);
            }
            f /= 2;
        }
    }
};

#endif
