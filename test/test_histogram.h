#ifndef test_histogram_h
#define test_histogram_h

#include "gtest/gtest.h"
#include "histogram.h"


TEST(Histogram, int) {
    histogram::Histogram<int> histogram(0, 1024, 1024);

    for (unsigned int i = 0; i <= 1024; i++)
        ASSERT_EQ(i, histogram.bin(i));
    ASSERT_EQ(0, histogram.bin(-1));
    ASSERT_EQ(1024, histogram.bin(1025));

    for (unsigned int i = 0; i <= 1024; i++)
        ASSERT_EQ(i, histogram.value(i));

    for (unsigned int i = 1; i <= 1024; i++) {
        ASSERT_EQ(i, histogram.value(histogram.bin(i)));
        ASSERT_EQ(i - 1, histogram.value(histogram.bin(i) - 1));
    }

    histogram.add(1025);
    ASSERT_EQ(histogram.bin(1025), 1024);
    ASSERT_EQ(histogram[1024], 1);
}


TEST(Histogram, float) {
    histogram::Histogram<double> histogram(0, 1024, 1024);

    for (unsigned int i = 0; i <= 1024; i++)
            ASSERT_EQ(i, histogram.bin(i));
    ASSERT_EQ(0, histogram.bin(-1));
    ASSERT_EQ(1024, histogram.bin(1025));

    for (unsigned int i = 0; i <= 1024; i++)
            ASSERT_EQ(i, histogram.value(i));

    for (unsigned int i = 0; i < 1024; i++) {
        ASSERT_EQ(i, histogram.value(histogram.bin(i)));
    }
}

TEST(Histogram, floatUniform) {
    histogram::Histogram<double> histogram(0, 1, 10);

    ASSERT_EQ(0, histogram.bin(-0.1));
    ASSERT_EQ(0, histogram.bin(0));
    ASSERT_EQ(0, histogram.bin(0.0000001));
    ASSERT_EQ(0, histogram.bin(0.0999999));
    ASSERT_EQ(5, histogram.bin(0.5000001));
    ASSERT_EQ(5, histogram.bin(0.5999999));
    ASSERT_EQ(9, histogram.bin(0.9000001));
    ASSERT_EQ(9, histogram.bin(0.9999999));
    ASSERT_EQ(10, histogram.bin(1));
    ASSERT_EQ(10, histogram.bin(1.1));

    EXPECT_NEAR(0.1, histogram.value(histogram.bin(0.1000001)), 0.00001);
    EXPECT_NEAR(0.1, histogram.value(histogram.bin(0.1999999)), 0.00001);
    EXPECT_NEAR(0.5, histogram.value(histogram.bin(0.5000001)), 0.00001);
    EXPECT_NEAR(0.5, histogram.value(histogram.bin(0.5999999)), 0.00001);
}


TEST(Log2Histogram, int) {
    histogram::Log2Histogram<int> histogram(1, 1024, log2(1024));

    for (unsigned int i = 0; i <= 10; i++)
            ASSERT_EQ(i, histogram.bin((int)pow(2, i)));
    ASSERT_EQ(0, histogram.bin(-1));
    ASSERT_EQ(10, histogram.bin(1025));

    for (unsigned int i = 0; i <= 10; i++)
            ASSERT_EQ(pow(2, i), histogram.value(i));

    for (unsigned int i = 0; i <= 10; i++) {
        ASSERT_EQ((int)pow(2, i), histogram.value(histogram.bin((int)pow(2, i))));
    }
}


TEST(FloatLog2Histogram, float) {
    histogram::Log2Histogram<double> histogram(1, 1024, log2(1024));

    for (unsigned int i = 0; i <= 10; i++)
        ASSERT_EQ(i, histogram.bin((int)pow(2, i)));
    ASSERT_EQ(0, histogram.bin(-1));
    ASSERT_EQ(10, histogram.bin(1025));

    for (unsigned int i = 0; i <= 10; i++)
        ASSERT_EQ(pow(2, i), histogram.value(i));

    for (unsigned int i = 0; i <= 10; i++)
        ASSERT_EQ(pow(2, i), histogram.value(histogram.bin((int)pow(2, i))));
}

#endif
