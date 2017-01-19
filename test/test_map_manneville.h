#ifndef test_map_manneville_h
#define test_map_manneville_h

#include "map.h"
#include "gtest/gtest.h"

TEST(TestManneville, iteration) {
    map::Manneville map;

    Vector point(1);
    point[0] = "0.5";
    Vector expected_result(1);
    expected_result[0] = "0.75";
    map.T(point);
	EXPECT_EQ(point, expected_result);
}

TEST(TestManneville, derivative) {
    map::Manneville map;

    Vector point(1);
    point[0] = "0.5";
    Vector vector(1);
    vector[0] = "1";

    Vector expected_result(1);
    expected_result[0] = "2";

    map.dT(point, vector);

	EXPECT_EQ(expected_result[0], vector[0]);
}

#endif
