/*
Compares the expected divergence (from Lyapunov) with the measured divergence (from distance)
*/
#include "map.h"
#include "sampler.h"


// a observable that also stores the final state
class LyapunovFinal : public observable::Lyapunov {
public:
    Vector final_state;
    LyapunovFinal(map::Map & map, unsigned int tobs) : Lyapunov(map, tobs) {}

    virtual void finalise(Vector const& point) {
        observable::Lyapunov::finalise(point);
        final_state = point;
    }

    LyapunovFinal & operator=(LyapunovFinal const& other) {
        observable::Lyapunov::operator=(other);
        this->final_state = other.final_state;
        return *this;
    }
};


class TestHistogram : public SamplingHistogram<LyapunovFinal> {
    typedef LyapunovFinal Observable;
protected:
    std::vector<double> log_d_histogram;
    std::vector<double> log_d2_histogram;
    std::vector<std::vector<double> > all;
    Observable observable;
    Float sigma0;

    void export_log_d(std::string file_name, std::string directory="") const {

        std::vector<std::vector<double> > data;
        for (unsigned int bin = 0; bin <= this->bins(); bin++) {
            unsigned int count = (*this)[bin];
            if (count > 0) {
                std::vector<double> row(3);
                row[0] = this->value(bin);
                row[1] = log_d_histogram[bin]*1./count;
                // variance
                row[2] = log_d2_histogram[bin]*1./count;
                row[2] -= row[1]*row[1];
                // std
                row[2] = sqrt(row[2]);

                data.push_back(row);
            }
        }

        io::save(data, directory + file_name);
    }

public:

    TestHistogram(double lowerBound, double upperBound, unsigned int bins, Observable const& observable, Float sigma0)
            : SamplingHistogram<LyapunovFinal>(lowerBound, upperBound, bins),
              log_d_histogram(this->bins()),
              log_d2_histogram(this->bins()),
              observable(observable), sigma0(sigma0)
    {}

    virtual void measure(Observable const& result, Observable const& result_prime, double acceptance) {
        SamplingHistogram<LyapunovFinal>::measure(result, result_prime, acceptance);

        Float sigma = sigma0/result.stretch();
        Observable result2(observable);
        result2.observe(result.state + sigma*result.eigenvector());

        Float dist = (result2.final_state - result.final_state).norm();
        double value = (dist/sigma0).toDouble();

        log_d_histogram[this->bin(result.observable())] += value;
        log_d2_histogram[this->bin(result.observable())] += value*value;

        std::vector<double> row(3);
        row[0] = result.lyapunov();
        row[1] = result2.lyapunov();
        row[2] = value;
        all.push_back(row);
    }

    virtual void export_histogram(std::string file_name, std::string directory="") const {
        SamplingHistogram<LyapunovFinal>::export_pretty(file_name, directory);

        export_log_d("log_d_" + file_name, directory);

        io::save(all, directory + "all_" + file_name);
    }

    void reset() {
        SamplingHistogram<LyapunovFinal>::reset();

        std::fill(log_d_histogram.begin(), log_d_histogram.end(), 0);
        std::fill(log_d2_histogram.begin(), log_d2_histogram.end(), 0);
    }
};


void tent3_us_up(unsigned int tobs, unsigned int samples, int log_sigma0, std::string directory="./") {
    mpfr::mpreal::set_default_prec(256);

    map::Tent map(3);
    LyapunovFinal observable(map, tobs);

    Float sigma0 = exp(-Float(log_sigma0));
    std::cout << "sigma0 = " << sigma0 << std::endl;

    TestHistogram histogram(log(3./2) - 0.0001, log(3) - 0.0001, tobs, observable, sigma0);
    proposal::Uniform<LyapunovFinal> proposal(map.boundary);

    MetropolisHastings<LyapunovFinal> mc(observable, proposal, histogram);

    mc.sample(samples);

    histogram.export_histogram(format("divergence_tent3_us_up_%d_%d.dat", tobs, log_sigma0), directory);
}


void sm6_us_up(unsigned int tobs, unsigned int samples, int log_sigma0, std::string directory="./") {
    mpfr::mpreal::set_default_prec(256);

    map::Standard map(6);
    LyapunovFinal observable(map, tobs);

    Float sigma0 = exp(-Float(log_sigma0));
    std::cout << "sigma0 = " << sigma0 << std::endl;

    TestHistogram histogram(0, log(4) - 0.00001, tobs, observable, sigma0);
    proposal::Uniform<LyapunovFinal> proposal(map.boundary);

    MetropolisHastings<LyapunovFinal> mc(observable, proposal, histogram);

    mc.sample(samples);

    histogram.export_histogram(format("divergence_sm6_us_up_%d_%d.dat", tobs, log_sigma0), directory);
}


void sm6_wl_sy(unsigned int tobs, unsigned int samples, int log_sigma0, std::string directory="./") {
    mpfr::mpreal::set_default_prec(256);

    map::Standard map(6);
    LyapunovFinal observable(map, tobs);

    Float sigma0 = exp(-Float(log_sigma0));
    std::cout << "sigma0 = " << sigma0 << std::endl;
    TestHistogram histogram(0, log(4) - 0.00001, tobs, observable, sigma0);
    proposal::PowerLawIsotropic<LyapunovFinal> proposal(map.boundary, -3, 40);

    WangLandau<LyapunovFinal> mc(observable, proposal, histogram);

    mc.sample(10, samples);

    histogram.export_histogram(format("divergence_sm6_wl_sy_%d_%d.dat", tobs, log_sigma0), directory);
}


void tent3_wl_sy(unsigned int tobs, unsigned int samples, int log_sigma0, std::string directory="./") {
    mpfr::mpreal::set_default_prec(256);

    map::Tent map(3);
    LyapunovFinal observable(map, tobs);

    Float sigma0 = exp(-Float(log_sigma0));
    std::cout << "sigma0 = " << sigma0 << std::endl;
    TestHistogram histogram(0, log(4) - 0.00001, tobs, observable, sigma0);
    proposal::PowerLawIsotropic<LyapunovFinal> proposal(map.boundary, -3, 40);

    WangLandau<LyapunovFinal> mc(observable, proposal, histogram);

    mc.sample(10, samples);

    histogram.export_histogram(format("divergence_tent3_wl_sy_%d_%d.dat", tobs, log_sigma0), directory);
}


int main() {
    // varying t and fixed \Delta
    tent3_us_up(10, 10000, 5);
    tent3_us_up(20, 10000, 5);
    tent3_us_up(40, 10000, 5);

    sm6_us_up(10, 10000, 5);
    sm6_us_up(20, 10000, 5);
    sm6_us_up(40, 10000, 5);

    // fixed t and varying \Delta
//    sm6_us_up(20, 10000, 3);
//    sm6_us_up(20, 10000, 4);
//    sm6_us_up(20, 10000, 5);
//    sm6_us_up(20, 10000, 6);
//    sm6_us_up(20, 10000, 7);
    return 0;
}
