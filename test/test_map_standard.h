#ifndef test_map_standard_h
#define test_map_standard_h

#include "map.h"

TEST(TestStandardMap, iteration) {
    map::Standard map(Float(4));

    Vector point(2);
    point << Float("0.1"), Float("0.1");

    map.T(point);

    // x = 1/10 + 4/(2*pi)*sin(2*pi/10); y = 1/10 + x
    Vector expected_result(2);
    expected_result << Float("0.47419571351545561328516762646"),
                       Float("0.57419571351545561328516762646");

    for (unsigned int i = 0; i < expected_result.size(); i++) {
        EXPECT_DOUBLE_EQ((double)expected_result[i], (double)point[i]);
    }
}

TEST(TestStandardMap, derivative) {
    map::Standard map(Float(4));

    Vector point(2);
    point << Float("0.1"), Float("0.1");
    Vector vector(2);
    vector << 1/sqrt(Float(2)), 1/sqrt(Float(2));

    map.dT(point, vector);

    // x = 1/10 + 4/(2*pi)*sin(2*pi/10); y = 1/10 + x
    Vector expected_result(2);
    expected_result << Float("2.995352392457285"), Float("3.702459173643832");

    for (unsigned int i = 0; i < expected_result.size(); i++) {
        EXPECT_DOUBLE_EQ((double)expected_result[i], (double)vector[i]);
    }
}


#endif
