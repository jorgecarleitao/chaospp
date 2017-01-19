# Chaos++

Chaos++ is a headers-only C++11 library to sample chaotic systems. 
Its main goal is to provide a common library to apply importance sampling Monte Carlo techniques 
to study chaotic systems.

Main publication:

* (1) To be published.

Relevant references:

* (2) [Monte Carlo Sampling in Fractal Landscapes](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.220601)
* (3) [Efficiency of Monte Carlo sampling in chaotic systems](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.052916)
* (4) [Stagger-and-Step Method: Detecting and Computing Chaotic Saddles in Higher Dimensions](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.2261)
* (5) [Large deviations of Lyapunov exponents](http://iopscience.iop.org/article/10.1088/1751-8113/46/25/254002/meta)
* (6) [Precision shooting: Sampling long transition pathways](http://scitation.aip.org/content/aip/journal/jcp/129/19/10.1063/1.2978000)
* (7) [Multicanonical MCMC for sampling rare events: an illustrative review](http://link.springer.com/article/10.1007%2Fs10463-014-0460-2)

## Dependencies

Chaos++ depends on gmp, mpfr, and Eigen. gmp+mpfr are libraries that allow to use arbitrary precision.
Eigen is a library for algebra operations, which Chaos++ uses for Jacobian matrixes, eigenvalues, and QR decomposition.

## Installation

Chaos++ is a header-only library. You only need to add "chaospp" directory to your include path and `#include` the relevant 
headers in your c++ code.

To compile the code, `include` the dependencies and add the library paths to them. This package contains a
[cmake](https://cmake.org/) to facilitate the compilation process across different operating systems.

## Examples


### Sampling 
Compute the distribution of FTLE (FT=10) of the Standard map with K=6:

    #include "map.h"
    #include "sampler.h"

    // define the map (K=6, SM6)
    map::Standard map(6);

    // define the observable (FTLE with FT=10)
    typedef observable::Lyapunov Obs;
    Obs observable(map, 10);

    // define histogram boundaries and number of bins
    // from [0, log(6)] with 100 bins
    SamplingHistogram<Obs> histogram(0 - 0.00001, log(6) - 0.00001, 100);

    // select proposal (uniform proposal)
    proposal::Uniform<Obs> proposal(map.boundary);

    // Select sampler (Metropolis-Hastings or WangLandau)
    MetropolisHastings<Obs> mc(observable, proposal, histogram);

    // sample
    mc.sample(10000); // 10000 samples

    // export histogram
    histogram.export_pretty("s.dat");

To use another distribution, e.g. canonical ensemble, use:

    double beta = 1;
    for(unsigned int bin = 0; bin < histogram.log_pi.size(); bin++) {
        histogram.log_pi[bin] = -beta*bin;
    }

before calling `mc.sample`.

### Find states 

Find a state with a large escape time:

    #include "map.h"
    #include "optimizer.h"

    // define the map (K=6, SM6)
    map::Standard map(6);

    // define the optimizer, target escape time = 20
    optimizer::Adaptive optimizer(map, 20);

    // call it to get a point
    std::cout << optimizer.get_point().state << std::endl;

Numerous examples, which reproduce most of the results published in Refs. (1-3), are available in directories
`examples`, `sample`, `search`, `test_assumptions`, and `test_canonical`.
Each file `*.cpp` is an example that obtains what is described inside that file. 

## Tests

This package includes a set of unit tests in directory `test` to validate relevant parts of the code.
These are written in [googletest](https://github.com/google/googletest).

## Functionality implemented

### Maps

* Standard map
* Tent map
* Open tent map
* 2 coupled Standard maps
* Pommeau-Manneville map
* N-coupled Henon maps

Implemented evolution both in phase-space and in tangent space, in arbitrary precision.

(defined in `map.h`)

### Proposals

* Uniform proposal
* Power-Law Isotropic proposal
* Lyapunov Isotropic proposal
* tstar proposal
* Adaptive proposal
* Anisotropic proposal

(defined in `proposal.h`)

### Algorithms/Methodologies

* Metropolis-Hastings algorithm (arbitrary target distribution)
* Wang-Landau algorithm (converges to MH with flat-histogram)
* Hill climbing (maximize/minimize)

(defined in `sampling.h` and `optimization.h`)

### Observables

* Escape time (`observable::EscapeTime`)
* FT Lyapunov exponent (`observable::Lyapunov`)

(defined in `observables.h`)

For example,

    map::OpenTent map(3, 5);
    observable::EscapeTime obs(map);

    Vector point(1);
    point[0] = Float("0.0000000001");

    obs.observe(point);

    // `obs.escape_time` contains the escape time of the state.


### Other functionality

* an histogram template class to create histograms to both discrete and continuous variables (`histogram.h`)
* functions to import and export arbitrary std::vector's as TSV or CSV (`io.h`)
* other math funtionality.

(defined in `histogram.h`, `io.h` and `auxiliar.h`)

## FAQ

### How to measure other quantities at each markov step?

The class `SamplingHistogram` has a method `measure` that is called by `MetropolisHastings` after each markov step. 
You can subclass `SamplingHistogram` and overload `measure`. For example, to measure and export the acceptance,
use:

    // An histogram for escape time that also measures the acceptance rate
    class HistogramWithAcceptance : public SamplingHistogram<observable::EscapeTime> {
    protected:
        std::vector<double> acceptance_histogram;

        void export_acceptance(std::string file_name, std::string directory="") const {
    
            std::vector<std::vector<double> > data;
            for (unsigned int bin = 0; bin <= this->bins(); bin++) {
                unsigned int count = (*this)[bin];  // count of the bin
                if (count > 0) {
                    std::vector<double> row(2);
                    row[0] = this->value(bin);  // the value of the bin
                    row[1] = acceptance_histogram[bin]*1./count;  // mean acceptance
                    data.push_back(row);
                }
            }
    
            io::save(data, directory + file_name);
        }
    
    public:

        HistogramWithAcceptance(unsigned int lowerBound, unsigned int upperBound, unsigned int bins)
                : SamplingHistogram<observable::EscapeTime>(lowerBound, upperBound, bins),
                  acceptance_histogram(this->bins()) {}
    
        void measure(observable::EscapeTime const& result, observable::EscapeTime const& result_prime, double acceptance) {
            // measure the escape time histogram, as before,
            SamplingHistogram<observable::EscapeTime>::measure(result, result_prime, acceptance);
            
            // And also store acceptance...
            acceptance_histogram[this->bin(result.observable())] += acceptance;
        }

        void export_histogram(std::string file_name, std::string directory="") const {
            SamplingHistogram<observable::EscapeTime>::export_pretty(file_name, directory);
            export_acceptance("acceptance_" + file_name, directory);
        }

        void reset() {
            SamplingHistogram<observable::EscapeTime>::reset();
    
            // reset to 0
            std::fill(acceptance_histogram.begin(), acceptance_histogram.end(), 0);
        }
    };

### How to measure other quantities of the state/at each map iteration?

The class `observable::Observable` has a method `observe(state)` that is called by `MetropolisHastings` to compute
anything about that state. For example, `observable::EscapeTime` evolves the state until it exits the system and stores
its escape time.
Subclass `observable::Observable` or any of its subclasses to make new measurements about the state. 
Subclass `observable::EscapeTime` and oberload `evolve` to make new measurements at each time-step.
In `observables.h` you can find examples where this is done, e.g. to also compute the evolution in the tangent space.
You then use the new observable in your code, using e.g. in the example above, use

    typedef MyObservable Obs;

## Authors

This code was written by [Jorge C. Leitão](http://jorgecarleitao.net).
