#ifndef chaospp_optimization_h
#define chaospp_optimization_h

#include <vector>

#include "proposal.h"
#include "map.h"
#include "observable.h"

namespace optimizer {

class Profiler {
    typedef observable::EscapeTime Observable;
public:
    virtual void start(Observable const& result) = 0;
    virtual void measure(Observable const& result, Observable const& result_prime, Float const& delta, double acceptance) = 0;

    virtual void export_all(std::string directory, std::string file_name) const = 0;
};


template <typename Observable>
class Optimizer {
private:
    //! a list of profilers that store information
    std::vector<Profiler*> profilers;

    void measure(Observable const& result, Observable const& result_prime, Float const& delta) {
        for (auto* profiler : profilers) {
            profiler->measure(result, result_prime, delta, 1);
        }
    }

    void start_profilers(Observable const& result) {
        for (auto* profiler : profilers) {
            profiler->start(result);
        }
    }
protected:
    unsigned int max_time;
    proposal::Proposal<Observable> & proposal;
    Observable observable;
public:

    Optimizer(Observable const& observable, proposal::Proposal<Observable> & proposal, unsigned int max_time) :
            observable(observable), proposal(proposal), max_time(max_time) {}

    virtual Observable get_point(unsigned int max_trials = 0) {
        Observable result(this->observable);
        result.observe(this->proposal.proposeUniform());
        start_profilers(result);

        unsigned int trial = 0;
        while (result.escape_time < max_time and (max_trials == 0 or trial < max_trials)) {
            trial++;
            Observable result_prime(this->observable);
            result_prime.observe(this->proposal.propose(result));

            this->measure(result, result_prime, proposal.get_delta());
            proposal.update(result, result_prime);
            if (result_prime.escape_time > result.escape_time)
                trial = 0;
            if (result_prime.escape_time >= result.escape_time) {
                result = result_prime;
            };
        }
        return result;
    }

    //! Adds a profiler to the list of profilers. Its `update()` will be called on every search.
    void add_profiler(Profiler & profiler) {
        profilers.push_back(&profiler);
    }
};


//! Optimization using Power-Law proposal.
class PowerLaw : public Optimizer<observable::EscapeTime> {
public:
    PowerLaw(map::Map & map, unsigned int max_time, double min_s, double max_s) :
        Optimizer<observable::EscapeTime>(observable::EscapeTime(map), *new proposal::PowerLawIsotropic<observable::EscapeTime>(map.boundary, min_s, max_s), max_time) {}
};


//! Optimization using adaptive proposal.
class Adaptive : public Optimizer<observable::EscapeTime> {
public:
    Adaptive(map::Map & map, unsigned int max_time) :
        Optimizer<observable::EscapeTime>(observable::EscapeTime(map, max_time), *new proposal::Adaptive<observable::EscapeTime>(map.boundary), max_time) {}
};


//! Optimization using Isotropic Lyapunov proposal
class Isotropic : public Optimizer<observable::EscapeWithVector> {
public:
    Isotropic(map::Map & map, unsigned int max_time) :
            Optimizer<observable::EscapeWithVector>(observable::EscapeWithVector(map, max_time), *new proposal::LyapunovIsotropic<observable::EscapeWithVector>(map.boundary), max_time) {}
};


//! Optimization using anisotropic proposal.
class Anisotropic : public Optimizer<observable::EscapeWithMatrix> {
public:
    Anisotropic(map::Map & map, unsigned int max_time) :
        Optimizer<observable::EscapeWithMatrix>(observable::EscapeWithMatrix(map, max_time), *new proposal::Anisotropic<observable::EscapeWithMatrix>(map.boundary), max_time) {}
};

}

#endif
