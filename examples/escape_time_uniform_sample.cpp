#include "map.h"
#include "sampler.h"

typedef observable::EscapeTime Obs;


void measure_us(map::Map & map, unsigned int samples) {
    Obs observable(map, 200);

    SamplingHistogram<Obs> histogram(0, 200, 200);
    proposal::Uniform<Obs> proposal(map.boundary);

    MetropolisHastings<Obs> mc(observable, proposal, histogram);

    mc.sample(samples);

    histogram.export_histogram(format("%s_us_%d.dat", map.name.c_str(), samples));
}


void measure_cs(map::Map & map, unsigned int samples) {
    Obs observable(map, 200);

    SamplingHistogram<Obs> histogram(0, 200, 200);
    proposal::Uniform<Obs> proposal(map.boundary);

    // set log_pi to be the canonical ensemble, proportional to -E*\beta
    double beta = -1;  // -1 favors larger values of E(=t_e here)
    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -beta*bin;
    }

    MetropolisHastings<Obs> mc(observable, proposal, histogram);

    mc.sample(samples, samples);

    histogram.export_histogram(format("%s_cs_%d.dat", map.name.c_str(), samples));
}


int main() {
    mpfr::mpreal::set_default_prec(64);

    map::NCoupledHenon ncoupled(4);
    measure_us(ncoupled, 100000);

    measure_cs(ncoupled, 100000);
    return 0;
}
