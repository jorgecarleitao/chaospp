#include "gtest/gtest.h"

#include "test_sampling.h"
#include "test_observables.h"
#include "test_map_manneville.h"
#include "test_map_standard.h"
#include "test_map_ncoupled.h"
#include "test_histogram.h"
#include "test_isotropic_proposal.h"
#include "test_anisotropic_proposal.h"


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
