#ifndef escape_test_ncoupled_h
#define escape_test_ncoupled_h

#include "observable.h"
#include "map.h"


TEST(NCoupledHenon, T) {
    map::NCoupledHenon map(6);

    Vector point(6);
    point[0] = "0.11";
    point[1] = "0.11";
    point[2] = "0.11";
    point[3] = "0.11";
    point[4] = "0.11";
    point[5] = "0.11";
    map.T(point);
    EXPECT_DOUBLE_EQ((double)Float("3.0209"), (double)point[0]);
    EXPECT_DOUBLE_EQ((double)Float("4.0209"), (double)point[1]);
    EXPECT_DOUBLE_EQ((double)Float("5.0209"), (double)point[2]);
}

TEST(NCoupledHenon, escape_time) {
    map::NCoupledHenon map(6);
    observable::EscapeTime observer(map);

    Vector point(6);
    point[0] = "0.11";
    point[1] = "2.15";
    point[2] = "2.11";
    point[3] = "0.11";
    point[4] = "0.11";
    point[5] = "0.11";

    observer.observe(point);

    ASSERT_EQ(3, observer.escape_time);
}

TEST(NCoupledHenon, jacobian2D) {
    map::NCoupledHenon map(2);
    observable::EscapeTime observer(map);

    Vector point(2);
    point[0] = "0.11";
    point[1] = "0.11";
    Matrix jac = map.jacobian(point);
    ASSERT_EQ(Float("-0.22"), jac(0, 0));
    ASSERT_EQ(Float("0.3"), jac(0, 1));
    ASSERT_EQ(Float("1"), jac(1, 0));
    ASSERT_EQ(Float("0"), jac(1, 1));
}

TEST(NCoupledHenon, jacobian4D) {
    map::NCoupledHenon map(4);

    Vector point(4);
    point[0] = "0.11";
    point[1] = "2.15";
    point[2] = "0.11";
    point[3] = "0.11";
    Matrix jac = map.jacobian(point);

    EXPECT_DOUBLE_EQ((double)Float("0.18"), (double)jac(0, 0)); // 2*x + k
    EXPECT_DOUBLE_EQ((double)Float("-0.4"), (double)jac(0, 1));  // B
    EXPECT_DOUBLE_EQ((double)Float("-0.4"), (double)jac(1, 0));  // B
    EXPECT_DOUBLE_EQ((double)Float("-3.9"), (double)jac(1, 1)); // 2*x + k
    EXPECT_DOUBLE_EQ((double)Float("0.3"), (double)jac(0, 2)); // -k
    EXPECT_DOUBLE_EQ((double)Float("0.3"), (double)jac(1, 3)); // -k
    EXPECT_DOUBLE_EQ((double)Float("1"), (double)jac(2, 0)); // 1
    EXPECT_DOUBLE_EQ((double)Float("1"), (double)jac(3, 1)); // 1
}


#endif
