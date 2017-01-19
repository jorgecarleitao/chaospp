#include "map.h"
#include "sampler.h"


void measure(map::Map & map, std::string map_name, double delta, std::string directory="./") {
    mpfr::mpreal::set_default_prec(64);

    proposal::LyapunovIsotropic<observable::EscapeWithVector> proposal(map.boundary, delta);

    observable::EscapeWithVector result(map);
    observable::EscapeWithVector result_prime(map);

    std::vector<std::vector<double> > data;
    for (unsigned int i = 0; i < 1000000; i++) {
        Vector state = proposal.proposeUniform();

        result.observe(state);

        Vector state_prime = proposal.propose(result);

        result_prime.observe(state_prime);

        std::vector<double> temp = {result.state[0].toDouble(), (double)result.escape_time, result.lyapunov(),
                                    result_prime.state[0].toDouble(), (double)result_prime.escape_time, result_prime.lyapunov(),
                                    proposal.log_acceptance(result, result_prime)};
        data.push_back(temp);
    }

    io::save(data, directory + map_name + format("_fh_%.2f", delta) + ".dat");
}


int main() {
    map::NCoupledHenon map(16);
    std::string map_name(format("ncoupled%d", 16));

    measure(map, map_name, 10);
    std::cout << 16 << std::endl;
}
