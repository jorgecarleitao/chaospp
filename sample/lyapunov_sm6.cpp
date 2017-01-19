/*
 Computes the finite-time Lyapunov exponent, acceptance and |x-x'| of the proposals
*/

#include <algorithm> // for max

#include "map.h"
#include "sampler.h"
#include "thesis_proposal.h"


typedef observable::Lyapunov Obs;


class TestHistogram : public SamplingHistogram<Obs> {
protected:
    std::vector<double> acceptance;
    std::vector<double> acceptance2;

    std::vector<double> delta_e;
    std::vector<double> delta_e2;

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

public:
    TestHistogram(double lowerBound, double upperBound, unsigned int bins)
            : SamplingHistogram<Obs>(lowerBound, upperBound, bins),
              acceptance(this->bins()),
              acceptance2(this->bins()),
              delta_e(this->bins()),
              delta_e2(this->bins()) {}

    virtual void measure(Obs const& result, Obs const& result_prime, double accept) {
        SamplingHistogram<Obs>::measure(result, result_prime, accept);
        unsigned int bin = this->bin(result.observable());

        acceptance[bin] += accept;
        acceptance2[bin] += accept*accept;

        double value = result.tobs*(result_prime.lyapunov() - result.lyapunov());
        delta_e[bin] += value;
        delta_e2[bin] += value*value;
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        SamplingHistogram<Obs>::export_pretty(file_name, directory);
        export_acceptance("acceptance_" + file_name, directory);
        export_delta_e("deltae_" + file_name, directory);
    }

    void reset() {
        SamplingHistogram<Obs>::reset();

        // reset to 0
        std::fill(acceptance.begin(), acceptance.end(), 0);
        std::fill(acceptance2.begin(), acceptance2.end(), 0);
        std::fill(delta_e.begin(), delta_e.end(), 0);
        std::fill(delta_e2.begin(), delta_e2.end(), 0);
    }
};


// Tstar Proposal
void wl_tp(unsigned int tobs, double delta, std::string directory="./") {

    mpfr::mpreal::set_default_prec(64);

    map::Standard map(6);
    Obs observable(map, tobs);

    TestHistogram histogram(0, 1.8, 20);

    TstarProposal<Obs> proposal(map.boundary, delta, tobs, histogram);

    WangLandau<Obs> mc(observable, proposal, histogram);

    mc.sample(8, 100000);

    histogram.export_histogram(format("sm6_%d_wl_tp_%.2f.dat", tobs, delta), directory);
    histogram.export_entropy(format("sm6_%d_wl_tp_%.2f.dat", tobs, delta), directory);
}


// Uniform Proposal
void wl_up(unsigned int tobs, std::string directory="./") {

    mpfr::mpreal::set_default_prec(64);

    map::Standard map(6);
    Obs observable(map, tobs);

    proposal::Uniform<Obs> proposal(map.boundary);

    TestHistogram histogram(-1, 6, 10*tobs);

    WangLandau<Obs> mc(observable, proposal, histogram);

    mc.sample(8, 10000);

    histogram.export_histogram(format("sm6_%d_wl_up.dat", tobs), directory);
}


// Thesis Proposal in uniform sampling
void us_tp(unsigned int tobs, double delta, std::string directory="./") {

    mpfr::mpreal::set_default_prec(64);

    map::Standard map(6);
    Obs observable(map, tobs);

    TestHistogram histogram(-1, 6, 30);

    TstarProposal<Obs> proposal(map.boundary, delta, tobs, histogram);

    MetropolisHastings<Obs> mc(observable, proposal, histogram);

    mc.sample(100000);

    histogram.export_histogram(format("sm6_%d_us_tp_%.2f.dat", tobs, delta), directory);
}


// Uniform sampling
void us_up(unsigned int tobs, std::string directory="./") {

    mpfr::mpreal::set_default_prec(64);

    map::Standard map(6);
    Obs observable(map, tobs);

    TestHistogram histogram(-1, 6, 10*tobs);
    proposal::Uniform<Obs> proposal(map.boundary);

    MetropolisHastings<Obs> mc(observable, proposal, histogram);

    mc.sample(10000);

    histogram.export_histogram(format("sm6_%d_us_up.dat", tobs), directory);
    histogram.export_entropy(format("sm6_%d_us_up.dat", tobs), directory);
}


int main() {
    wl_tp(10, 1);
    std::cout << "10 ended" << std::endl;
    wl_tp(20, 1);
    std::cout << "20 ended" << std::endl;
    wl_tp(30, 1);

    wl_up(10);
    std::cout << "10 ended" << std::endl;
    wl_up(20);
    std::cout << "20 ended" << std::endl;
    wl_up(30);

    us_up(10);
    std::cout << "10 ended" << std::endl;
    us_up(20);
    std::cout << "20 ended" << std::endl;
    us_up(30);

    us_tp(10, 1);
    std::cout << "10 ended" << std::endl;
    us_tp(20, 1);
    std::cout << "20 ended" << std::endl;
    us_tp(30, 1);

    return 0;
}
