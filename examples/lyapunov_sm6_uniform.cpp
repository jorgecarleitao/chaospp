/*
Computes distribution of FTLE for different tobs in the standard map, using uniform sampling.
*/
#include "map.h"
#include "sampler.h"


void measure(unsigned int tobs) {

    map::Standard map(6);
    observable::Lyapunov observable(map, tobs);

    SamplingHistogram<observable::Lyapunov> histogram(0 - 0.00001, log(6) - 0.00001, 100);
    proposal::Uniform<observable::Lyapunov> proposal(map.boundary);

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(1000);

    histogram.export_pretty(format("sm6_%d_uniform.dat", tobs), "./results");
}

int main() {
    mpfr::mpreal::set_default_prec(64);

    measure(10);
    measure(20);
    measure(40);

    return 0;
}
