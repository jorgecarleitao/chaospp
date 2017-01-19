#include "io.h"
#include "map.h"
#include "proposal.h"


void measure(map::Map & map) {
    std::string file_name = format("sm%s.dat", map.name.c_str());

    mpfr::mpreal::set_default_prec(64);

    proposal::Uniform<observable::Lyapunov> proposal(map.boundary);

    std::vector<std::vector<double> > data;
    // 100 states
    for (unsigned int i = 0; i < 100; i++) {
        Vector state = proposal.proposeUniform();
        // iterated for 100 times
        for (unsigned int t = 0; t < 100; t++) {
            map.T(state);

            // store state
            std::vector<double> row(map.D);
            for (unsigned int d = 0; d < map.D; d++)
                row[d] = state[d].toDouble();
            data.push_back(row);
        }
    }
    io::save(data, "./results/" + file_name);
}

int main() {
    map::Standard map("6");

    measure(map);
    return 0;
}
