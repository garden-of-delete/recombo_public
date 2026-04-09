#include <gtest/gtest.h>
#include <algorithm>
#include <pseudorandom.h>

static int randVals[] = {
    43374145,
    992298682,
    205945772,
    1734434687,
    700760031,
    873520869,
    1067034391,
    11073440,
    1427712438,
    1728157535
};

static int randValsNumber = 10;

TEST(PseudorandomTest, ProducesDeterministicSequenceFromSeed)
{
    pseudorandom r;
    r.sRandSimple(42);
    for (int i = 0; i < randValsNumber; i++)
    {
        int roll = r.RandSimple();
        ASSERT_EQ(randVals[i], roll);
    }
}

TEST(PseudorandomTest, IntegerCoversExpectedRange)
{
   int low, hi, x;
   low = hi = 5;
   pseudorandom r(42);
   for(int i = 0; i < 1000; i++)
   {
      x = r.rand_integer(1, 11);
      low = std::min(low, x);
      hi = std::max(hi, x);
   }
   ASSERT_EQ(1, low);
   ASSERT_EQ(10, hi);
}
