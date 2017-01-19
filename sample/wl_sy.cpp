/*
This code compares Wang-Landau simulations with different proposals
*/

#include "map.h"
#include "sampler.h"


class TestHistogram : public SamplingHistogram<observable::EscapeWithVector> {
    typedef observable::EscapeWithVector Observable;
protected:
    std::vector<double> acceptance_histogram;
    std::vector<double> log_delta_histogram;
    std::vector<double> log_delta2_histogram;
    std::vector<double> conditional_histogram;
    std::vector<double> lambda_histogram;
    std::vector<double> lambda2_histogram;

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

    TestHistogram(unsigned int lowerBound, unsigned int upperBound, unsigned int bins)
            : SamplingHistogram<observable::EscapeWithVector>(lowerBound, upperBound, bins),
              acceptance_histogram(this->bins()),
              log_delta_histogram(this->bins()),
              log_delta2_histogram(this->bins()),
              conditional_histogram(this->bins()),
              lambda_histogram(this->bins()),
              lambda2_histogram(this->bins())
    {}

    virtual void measure(Observable const& result, Observable const& result_prime, double acceptance) {
        SamplingHistogram<observable::EscapeWithVector>::measure(result, result_prime, acceptance);
        unsigned int bin = this->bin(result.observable());

        acceptance_histogram[bin] += acceptance;

        if (0.1 < acceptance and acceptance < 0.9) {
            Float norm = aux::get_norm(result_prime.state - result.state);
            conditional_histogram[bin] += 1;
            double value = -log(norm).toDouble();
            log_delta_histogram[bin] += value;
            log_delta2_histogram[bin] += value * value;
        }
        double value = log(result.stretch()).toDouble();
        lambda_histogram[bin] += value;
        lambda2_histogram[bin] += value * value;
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        SamplingHistogram<observable::EscapeWithVector>::export_pretty(file_name, directory);
        export_acceptance("acceptance_" + file_name, directory);
        export_log_delta("log_delta_" + file_name, directory);
        export_conditional("conditional_acceptance_" + file_name, directory);
    }

    void reset() {
        SamplingHistogram<observable::EscapeWithVector>::reset();

        // reset to 0
        std::fill(acceptance_histogram.begin(), acceptance_histogram.end(), 0);
        std::fill(log_delta_histogram.begin(), log_delta_histogram.end(), 0);
        std::fill(log_delta2_histogram.begin(), log_delta2_histogram.end(), 0);
        std::fill(conditional_histogram.begin(), conditional_histogram.end(), 0);
        std::fill(lambda_histogram.begin(), lambda_histogram.end(), 0);
        std::fill(lambda2_histogram.begin(), lambda2_histogram.end(), 0);
    }
};

std::string directory = "./results/";


// wang_landau with Power-Law proposal
void standard_wl_sy(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::Standard map(6);

    TestHistogram histogram(0, 30, 30);
    proposal::PowerLawIsotropic<observable::EscapeWithVector> proposal(map.boundary, -2, 34);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("standard_wl_sy.dat", directory);
    histogram.export_entropy("standard_wl_sy.dat", directory);
}


// wang_landau with FTLE proposal
void standard_wl_lambda(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::Standard map(6);

    TestHistogram histogram(0, 30, 30);
    proposal::LyapunovIsotropic<observable::EscapeWithVector> proposal(map.boundary, 100);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("standard_wl_lambda.dat", directory);
    histogram.export_entropy("standard_wl_lambda.dat", directory);
}


// wang_landau with FTLE proposal
void ncoupled_wl_lambda(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::NCoupledHenon map(6);

    TestHistogram histogram(0, 30, 30);
    proposal::LyapunovIsotropic<observable::EscapeWithVector> proposal(map.boundary, 100);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("ncoupled_wl_lambda.dat", directory);
    histogram.export_entropy("ncoupled_wl_lambda.dat", directory);
}


// wang_landau with Power-Law proposal
void ncoupled_wl_sy(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::NCoupledHenon map(6);

    TestHistogram histogram(0, 30, 30);
    proposal::PowerLawIsotropic<observable::EscapeWithVector> proposal(map.boundary, -2, 34);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("ncoupled_wl_sy.dat", directory);
    histogram.export_entropy("ncoupled_wl_sy.dat", directory);
}


// wang_landau with FTLE proposal
void tent_wl_lambda(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::OpenTent map(3, 5);

    TestHistogram histogram(0, 30, 30);
    proposal::LyapunovIsotropic<observable::EscapeWithVector> proposal(map.boundary, 100);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("tent_wl_lambda.dat", directory);
    histogram.export_entropy("tent_wl_lambda.dat", directory);
}


// wang_landau with Power-Law proposal
void tent_wl_sy(unsigned int steps, unsigned int samples) {
    mpfr::mpreal::set_default_prec(128);

    map::OpenTent map(3, 5);

    TestHistogram histogram(0, 30, 30);
    proposal::PowerLawIsotropic<observable::EscapeWithVector> proposal(map.boundary, -2, 40);

    WangLandau<observable::EscapeWithVector> mc(map, proposal, histogram);

    mc.sample(steps, samples/steps);

    histogram.export_histogram("tent_wl_sy.dat", directory);
    histogram.export_entropy("tent_wl_sy.dat", directory);
}


int main() {
    //tent_wl_lambda(10, 40000);
    //standard_wl_lambda(10, 40000);
    //ncoupled_wl_lambda(10, 80000);

    tent_wl_sy(10, 400000);
    //standard_wl_sy(10, 40000);
    //ncoupled_wl_sy(10, 80000);
    return 0;
}
