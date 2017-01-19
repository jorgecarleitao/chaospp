#include "map.h"
#include "sampler.h"


int main() {
    char * env_max_t = std::getenv("MAX_ESCAPE_TIME");
    if (env_max_t == NULL) {std::cout << "MAX_ESCAPE_TIME not defined" << std::endl;exit(1);}
    unsigned int max_t = (unsigned int)atoi(env_max_t);

    char * env_samples = std::getenv("SAMPLES");
    if (env_samples == NULL) {std::cout << "SAMPLES not defined" << std::endl;exit(1);}
    unsigned int samples = (unsigned int)atoi(env_samples);

    char * env_D = std::getenv("DIMENSION");
    if (env_D == NULL) {std::cout << "DIMENSION not defined" << std::endl;exit(1);}
    unsigned int D = (unsigned int)atoi(env_D);

    mpfr::mpreal::set_default_prec(128);

    typedef SamplingHistogram<observable::EscapeWithVector> Histogram;
    typedef proposal::LyapunovIsotropic<observable::EscapeWithVector> Proposal;

    map::NCoupledHenon map(D);
    observable::EscapeWithVector observable(map);

    Histogram histogram(0, max_t, max_t);
    Proposal proposal(map.boundary, 10);

    WangLandau<observable::EscapeWithVector> mc(observable, proposal, histogram);

    observable.observe(proposal.proposeUniform());

    mc.approximate_entropy(10, samples);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < samples; i++) {
        histogram.reset();
        mc.round_trip(observable);

        std::vector<double> temp(1);
        temp[0] = histogram.count();
        data.push_back(temp);
    }

    io::save(data, format("results/round_trip_chm_d%d_t%d_s%d.dat",D, max_t, samples));
    return 0;
}
