#include <gtest/gtest.h>

#include <clkCigar.h>
#include <clkConformationBfacf3.h>
#include <clkConformationAsList.h>

using namespace std;

#define CRITICAL_Z (1.0 / 4.6852)  // Critical z-value for unknot from Schmirler (2012)

class BfacfRegressionTest : public ::testing::Test
{
protected:
    void TearDown() override
    {
        delete knot;
    }

    void initSingleComponent(int seed, double z)
    {
        clkCigar square;
        knot = new clkConformationBfacf3(square);
        knot->setSeed(seed);
        knot->setZ(z);
    }

    clkConformationBfacf3* knot = NULL;
};

TEST_F(BfacfRegressionTest, SingleComponentCriticalZ1000StepsSeed42)
{
    initSingleComponent(42, CRITICAL_Z);
    for (int i = 0; i < 1000; i++)
        knot->step();

    clkConformationAsList result(knot->getComponent(0));
    EXPECT_EQ(result.writeAsNewsud(), "dnus");
}

TEST_F(BfacfRegressionTest, SingleComponentCriticalZ5000StepsSeed2)
{
    initSingleComponent(2, CRITICAL_Z);
    for (int i = 0; i < 5000; i++)
        knot->step();

    clkConformationAsList result(knot->getComponent(0));
    EXPECT_EQ(result.writeAsNewsud(), "deusuwdn");
}

TEST_F(BfacfRegressionTest, SingleComponentCriticalZ1000StepsSeed42Reproducible)
{
    initSingleComponent(42, CRITICAL_Z);
    for (int i = 0; i < 1000; i++)
        knot->step();
    clkConformationAsList result1(knot->getComponent(0));
    string final1 = result1.writeAsNewsud();

    initSingleComponent(42, CRITICAL_Z);
    for (int i = 0; i < 1000; i++)
        knot->step();
    clkConformationAsList result2(knot->getComponent(0));
    string final2 = result2.writeAsNewsud();

    EXPECT_EQ(final1, final2);
}

TEST_F(BfacfRegressionTest, SingleComponentCriticalZ5000StepsSeed2Reproducible)
{
    initSingleComponent(2, CRITICAL_Z);
    for (int i = 0; i < 5000; i++)
        knot->step();
    clkConformationAsList result1(knot->getComponent(0));
    string final1 = result1.writeAsNewsud();

    initSingleComponent(2, CRITICAL_Z);
    for (int i = 0; i < 5000; i++)
        knot->step();
    clkConformationAsList result2(knot->getComponent(0));
    string final2 = result2.writeAsNewsud();

    EXPECT_EQ(final1, final2);
}

TEST_F(BfacfRegressionTest, DifferentSeedsDifferentTrajectories)
{
    initSingleComponent(1, CRITICAL_Z);
    for (int i = 0; i < 100; i++)
        knot->step();
    clkConformationAsList c1(knot->getComponent(0));
    string trajectory1 = c1.writeAsNewsud();

    initSingleComponent(2, CRITICAL_Z);
    for (int i = 0; i < 100; i++)
        knot->step();
    clkConformationAsList c2(knot->getComponent(0));
    string trajectory2 = c2.writeAsNewsud();

    EXPECT_NE(trajectory1, trajectory2);
}
