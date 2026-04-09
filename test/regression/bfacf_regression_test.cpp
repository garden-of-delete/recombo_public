#include <gtest/gtest.h>

#include <clkCigar.h>
#include <clkConformationBfacf3.h>

TEST(BfacfRegressionTest, Placeholder) {
    clkCigar square;
    clkConformationBfacf3 knot(square);
    EXPECT_EQ(knot.size(), 1);
}
