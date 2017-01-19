#include <iostream>

#include "map.h"
#include "optimizer.h"


int main() {
    mpfr::mpreal::set_default_prec(64);

    map::NCoupledHenon map(2);

    optimizer::Anisotropic optimizer(map, 20);

    std::cout << optimizer.get_point().state << std::endl;

    return 0;
}
