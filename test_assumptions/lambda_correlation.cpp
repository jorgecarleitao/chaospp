/*
 * Moves the trajectory t/2 steps and measures \lambda1, advances more t/2 and measures \lambda2,
 * output lambda1, lambda2. They should be independent of each other.
 */

#include "map.h"
#include "sampler.h"


void measure(std::string directory="./") {
    mpfr::mpreal::set_default_prec(64);

    /*
    typedef Map::Standard Map;
    Map map(6);
    std::string map_name("sm6");
    */
    map::Tent map(3);
    std::string map_name("tent3");
    /*
    typedef Map::Logistic Map;
    Map map(4);
    std::string map_name("logistic4");
    */

    proposal::Uniform<observable::EscapeTime> proposal(map.boundary);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < 10000; i++) {
        Vector state = proposal.proposeUniform();

        // measures [0, t/2]
        observable::Lyapunov result(map, 8);
        result.observe(state);

        // advances `state` to t/2
        for (unsigned int t = 0; t < 8; t++) {
            map.T(state);
        }

        // measures [t/2, t]
        observable::Lyapunov result_prime(map, 8);
        result_prime.observe(state);

        // stores both results
        std::vector<double> temp = {result.lyapunov(), result_prime.lyapunov()};
        data.push_back(temp);
    }

    // export results
    io::save(data, directory + map_name + ".dat");
}

int main() {
    measure();
}
