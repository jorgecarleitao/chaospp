/*
  Tests that a perturbation of \delta=e^-(lambda*t) guarantees that on average E(x+\delta) = E(x).
*/

#include "map.h"
#include "sampler.h"
#include "observable.h"


// Tests that lambda(x+\delta) = lambda(x).
void measure_lambda(std::string directory="./") {
    mpfr::mpreal::set_default_prec(64);

    /*
    typedef map::Standard Map;
    Map map(6);
    std::string map_name("sm6");
    */
    map::Tent map(3);
    std::string map_name("tent3");
    /*
    typedef map::Logistic Map;
    Map map(4);
    std::string map_name("logistic4");
    */
    std::string delta = "0.1";

    const static unsigned int tobs = 20;

    proposal::LyapunovIsotropic<observable::Lyapunov> proposal(map.boundary, delta);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < 100000; i++) {
        Vector state = proposal.proposeUniform();

        observable::Lyapunov result(map, tobs);
        result.observe(state);

        Vector state_prime = proposal.propose(result);

        observable::Lyapunov result_prime(map, tobs);
        result_prime.observe(state_prime);

        // stores both results
        std::vector<double> temp = {result.lyapunov(), result_prime.lyapunov() - result.lyapunov()};
        data.push_back(temp);
    }

    io::save(data, directory + format("lambda_%d_%s_%s.dat", tobs, map_name.c_str(), delta.c_str()));
}


// Tests that t_e(x+\delta) = t_e(x).
void measure_escape(std::string delta, std::string directory="./") {
    mpfr::mpreal::set_default_prec(64);
    std::cout << delta << std::endl;

    typedef map::Standard Map;
    Map map(6);
    std::string map_name("sm6");

//    typedef map::OpenTent Map;
//    Map map(3, 4);
//    std::string map_name("tent34");

    typedef observable::EscapeWithMatrix Observable;
    typedef proposal::LyapunovIsotropic<Observable> Proposal;

    Proposal proposal(map.boundary, delta);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < 10000; i++) {
        Vector state = proposal.proposeUniform();

        Observable result(map, 1024);
        result.observe(state);

        Vector state_prime = proposal.propose(result);

        Observable result_prime(map, 1024 + 24);
        result_prime.observe(state_prime);

        // stores both results
        std::vector<double> temp = {(double)result.escape_time,
                                    ((double) result_prime.escape_time) - ((double) result.escape_time)};
        data.push_back(temp);
    }

    io::save(data, directory + "escape_" + map_name + "_" + delta + ".dat");
}


int main() {
    measure_escape("1");
    measure_escape("0.1");
    measure_escape("0.01");
}
