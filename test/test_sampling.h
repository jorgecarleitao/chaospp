#ifndef chaospp_test_sampling_h
#define chaospp_test_sampling_h

#include "map.h"
#include "sampler.h"
#include "observable.h"

TEST(TestUniform, lyapunov_tent_map) {
    mpfr::mpreal::set_default_prec(64);

    map::Tent map(3);
    observable::Lyapunov observable(map, 10);

    SamplingHistogram<observable::Lyapunov> histogram(-2, 2, 400);
    proposal::Uniform<observable::Lyapunov> proposal(map.boundary);

    MetropolisHastings<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(1);

    EXPECT_EQ(histogram.count(), 1);
}


class SamplingHist : public SamplingHistogram<observable::EscapeTime> {
public:
    double mean_escape;

    SamplingHist(unsigned int lowerBound, unsigned int upperBound, unsigned int bins) : SamplingHistogram<observable::EscapeTime>(lowerBound, upperBound, bins) {
        mean_escape = 0;
    }

    virtual void measure(observable::EscapeTime const& result, observable::EscapeTime const& resultPrime, double acceptance) {
        SamplingHistogram<observable::EscapeTime>::measure(result, resultPrime, acceptance);
        // <x>_n += (x_n - <x>_{n-1})/n, see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
        mean_escape += (result.escape_time - mean_escape*1.)/count();
        //std::cout << result.escape_time << std::endl;
    }

    void reset() {
        SamplingHistogram<observable::EscapeTime>::reset();
        mean_escape = 0;
    }
};


// tests that we obtain the expected mean escape time under uniform sampling
TEST(TestUniform, escape_time_tent_map) {
    mpfr::mpreal::set_default_prec(64);

    map::OpenTent map(3, 5);
    observable::EscapeTime observable(map, 20);

    SamplingHist histogram(0, 20, 20);
    proposal::Uniform<observable::EscapeTime> proposal(map.boundary);

    MetropolisHastings<observable::EscapeTime> mc(observable, proposal, histogram);

    mc.sample(100000);

    double exp = 1./(1 - (1/3. + 1/5.));

    EXPECT_NEAR(histogram.mean_escape, exp, exp*0.05);
}


TEST(WangLandau, lyapunov_tent_map) {
    mpfr::mpreal::set_default_prec(64);

    map::Tent map(3);
    observable::Lyapunov observable(map, 10);

    SamplingHistogram<observable::Lyapunov> histogram(-2, 2, 400);
    proposal::PowerLawIsotropic<observable::Lyapunov> proposal(map.boundary, -1, 20);

    WangLandau<observable::Lyapunov> mc(observable, proposal, histogram);

    mc.sample(1, 50);

    EXPECT_EQ(histogram.count(), 50);
}


// tests that we obtain a flat-histogram using Wang-Landau and power-law distribution.
TEST(WangLandau, escape_time_tent_map) {
    mpfr::mpreal::set_default_prec(64);

    map::OpenTent map(3, 5);
    observable::EscapeTime observable(map, 10);

    SamplingHist histogram(0, 10, 10);
    proposal::PowerLawIsotropic<observable::EscapeTime> proposal(map.boundary, -1, 20);

    WangLandau<observable::EscapeTime> mc(observable, proposal, histogram);

    mc.sample(10, 10000);

    // The histogram is flat, thus the expected measured mean escape time is the mean
    // of the uniform distribution over [0, 10]
    double expected = (10 - 0)/2.;
    EXPECT_NEAR(histogram.mean_escape, expected, expected*0.01);
}


#endif
