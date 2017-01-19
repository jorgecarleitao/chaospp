#ifndef chaospp_test_anisotropic_proposal_h
#define chaospp_test_anisotropic_proposal_h

#include "gtest/gtest.h"
#include "proposal.h"
#include "observable.h"

struct EscapeWithTrajectory : public observable::EscapeWithMatrix {
    std::vector<Vector> trajectory;

    EscapeWithTrajectory(map::Map & map) : EscapeWithMatrix(map) {}

    void evolve(Vector & point) {
        EscapeWithMatrix::evolve(point);
        trajectory.push_back(point);
    }
};

TEST(Anisotropic, basic) {
    // Asserts that the anisotropic method generates an isotropic distribution at t = t_escape - 1.
    // Let \delta = r' - r, where r is a point with escape time 27 and r' is the proposal.
    // This test ensures that AVG[log10(|\delta[0]/\delta[1]|)] < 1, i.e. that the distance \delta on average is
    // close to isotropic.

    mpfr::mpreal::set_default_prec(512);

    map::NCoupledHenon map(4);
    EscapeWithTrajectory result(map);
    EscapeWithTrajectory resultPrime(map);

    Vector state(map.D);
    state << "2.247351146173699675939584601441840793919192562373547163157017081402955227531492710113525390625000000",
             "-1.141970787318847434230253208926983918540476764992444656116044043869806046131998300552368164062500000",
             "3.803983448890944066215028665089462570065497874736439552767475191785706556402146816253662109375000000",
             "1.083416859646563025778245376589602637749136687694059139053237572625221218913793563842773437500000000";
    result.observe(state);
    ASSERT_EQ(27, result.escape_time);

    proposal::Anisotropic<EscapeWithTrajectory> proposal(map.boundary);

    double avg = 0;
    for (unsigned int i = 0; i < 100; i++) {
        resultPrime.observe(proposal.propose(result));

        Vector end = result.trajectory[result.trajectory.size() - 1 - 1];
        Vector endPrime = resultPrime.trajectory[resultPrime.trajectory.size() - 1 - 1];
        Vector delta = endPrime - end;

        avg += log10(abs(delta[0]/delta[1])).toDouble();
    }
    avg /= 100;
    // assert the average is in the interval [-1, 1].
    ASSERT_NEAR(avg, 0, 1);
}

#endif
