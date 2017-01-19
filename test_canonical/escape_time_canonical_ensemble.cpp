#include "map.h"
#include "sampler.h"


class TestHistogram : public SamplingHistogram<observable::EscapeTime> {
protected:

    void export_acceptance(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (count > 0) {
                std::vector<double> row(2);
                row[0] = this->value(bin);
                row[1] = acceptance_histogram[bin]*1./count;
                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

public:
    std::vector<double> acceptance_histogram;

    TestHistogram(unsigned int lowerBound, unsigned int upperBound, unsigned int bins)
            : SamplingHistogram<observable::EscapeTime>(lowerBound, upperBound, bins),
              acceptance_histogram(this->bins()) {}

    virtual void measure(observable::EscapeTime const& result, observable::EscapeTime const& result_prime, double acceptance) {
        SamplingHistogram<observable::EscapeTime>::measure(result, result_prime, acceptance);
        acceptance_histogram[this->bin(result.observable())] += acceptance;
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        SamplingHistogram<observable::EscapeTime>::export_pretty(file_name, directory);
        export_acceptance("acceptance_" + file_name, directory);
    }

    void reset() {
        SamplingHistogram<observable::EscapeTime>::reset();

        // reset to 0
        std::fill(acceptance_histogram.begin(), acceptance_histogram.end(), 0);
    }
};


void canonical_uni(unsigned int samples, unsigned int thermalise_samples, double minus_beta, std::string directory="./") {
    mpfr::mpreal::set_default_prec(128);

    map::NCoupledHenon map(2);

    TestHistogram histogram(0, 40, 40);
    proposal::Uniform<observable::EscapeTime> proposal(map.boundary);

    // set log_pi to be the canonical ensemble, proportional to -E*\beta
    double beta = -minus_beta;  // -1 favors larger values of E(=t_e here)
    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -beta*bin;
    }

    MetropolisHastings<observable::EscapeTime> mc(map, proposal, histogram);

    mc.sample(samples, thermalise_samples);

    histogram.export_histogram(format("ce_uni_b%1.1f.dat", minus_beta), directory);
    histogram.export_entropy(format("ce_uni_b%1.1f.dat", minus_beta), directory);
}


void canonical_sy(unsigned int samples, unsigned int thermalise_samples, double minus_beta, std::string directory="./") {
    mpfr::mpreal::set_default_prec(128);

    map::NCoupledHenon map(2);

    TestHistogram histogram(0, 40, 40);
    proposal::PowerLawIsotropic<observable::EscapeTime> proposal(map.boundary, -3, 50);

    // set log_pi to be the canonical ensemble, proportional to -E*\beta
    double beta = -minus_beta;
    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -beta*bin;
    }

    MetropolisHastings<observable::EscapeTime> mc(map, proposal, histogram);

    mc.sample(samples, thermalise_samples);

    histogram.export_histogram(format("ce_sy_b%1.1f.dat", minus_beta), directory);
    histogram.export_entropy(format("ce_sy_b%1.1f.dat", minus_beta), directory);
}


int main() {

    // Uniform sampling
    canonical_uni(1000000, 0, 0);
    // canonical uniform proposal
    canonical_uni(900000, 10000, 1);
    // canonical Sweet-York proposal
    canonical_sy(900000, 10000, 1);
    return 0;
}
