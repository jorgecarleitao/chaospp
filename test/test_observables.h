#ifndef CHAOSPP_TEST_OBSERVABLES_H
#define CHAOSPP_TEST_OBSERVABLES_H

#include "observable.h"

TEST(TentMap, Lyapunov) {
    map::Tent map(3);
    observable::Lyapunov obs(map, 10);

    Vector point(1);
    point[0] = Float("0.0000000001");

    obs.observe(point);

    EXPECT_DOUBLE_EQ((double)log(Float(3)), obs.lyapunov());

    obs.observe(point);

    // check that copy copies correctly.
    observable::Lyapunov obs1(map, 10);
    obs1 = obs;

    EXPECT_DOUBLE_EQ(obs1.lyapunov(), obs.lyapunov());
}


TEST(TentMap, EscapeTime) {
    map::OpenTent map(3, 5);
    observable::EscapeTime obs(map);

    Vector point(1);

    // point in the exit set has escape time 1
    point[0] = Float("0.334");
    obs.observe(point);
    EXPECT_EQ(1, obs.escape_time);

    point[0] = Float("0.0000000001");
    obs.observe(point);
    EXPECT_EQ(21, obs.escape_time);

    // check that copy copies correctly.
    observable::EscapeTime obs1(map);
    obs1 = obs;

    EXPECT_EQ(21, obs1.escape_time);
}


#endif //CHAOSPP_TEST_OBSERVABLES_H
