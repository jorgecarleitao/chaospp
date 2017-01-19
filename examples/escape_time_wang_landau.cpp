/*
Example of the calculation of the distribution of escape time using Wang-Landau with isotropic proposal
*/

#include "map.h"
#include "sampler.h"

int main() {
    mpfr::mpreal::set_default_prec(64);

    map::NCoupledHenon map(4);
    observable::EscapeWithVector observable(map, 20);  // max_time 20

    SamplingHistogram<observable::EscapeWithVector> histogram(0, 20, 20);
    proposal::LyapunovIsotropic<observable::EscapeWithVector> proposal(map.boundary);

    WangLandau<observable::EscapeWithVector> mc(observable, proposal, histogram);

    mc.sample(10, 10000);

    histogram.export_histogram(format("%s_wl.dat", map.name.c_str()));
    histogram.export_entropy(format("%s_wl.dat", map.name.c_str()));

    return 0;
}
