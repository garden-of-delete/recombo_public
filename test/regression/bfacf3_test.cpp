// Extracted from test/unit/clk_test.cpp's Bfacf3TestFixture (2026-07-12).
// RandomReset and Bfacf3Run pin down the result of actually running the
// BFACF simulation (RNG determinism across reset/save/restore, and a
// literal target conformation after 100 steps) -- there is no closed-form
// answer to check these against, only "this is what a trusted version of
// the code produced." That makes them regression tests, not unit tests,
// which is why they moved here rather than staying in test/unit/.
//
// Bfacf3TestFixture's other two tests (Bfacf3, Bfacf3SetZ) stayed in
// clk_test.cpp -- they check behavior against independently-computable
// values (geometry, closed-form move probabilities), which is genuinely
// unit-shaped.
//
// Neither test here needs shared fixture state, so both are plain TEST()
// with a locally-constructed clkCigar rather than a fixture class.

#include <gtest/gtest.h>

#include <list>

#include <clkCigar.h>
#include <clkConformationAsList.h>
#include <clkConformationBfacf3.h>
#include <pseudorandom.h>
#include <threevector.h>

using namespace std;

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

TEST(Bfacf3RegressionTest, RandomReset)
{
   clkCigar square;
   int expected[] = {2, 46, 10, 80, 33, 41, 50, 1, 66, 80, 20, 71, 4, 41, 73, 71, 99, 70, 25, 4};
   list<int> l;

   // First verify that random number generator used by legacy is same as pseudorandom class.
   clkConformationBfacf3 knot(square);
   knot.setSeed(42);
   pseudorandom r(42);
   for (int i = 0; i < 20; i++)
   {
      int m = rand_integer(1, 100);
      int n = r.rand_integer(1, 100);
      ASSERT_EQ(expected[i], m);
      ASSERT_EQ(m, n);
      l.push_back(m);
   }

   // Next verify that resetting both random number generators with the same seed produces same result.
   r.sRandSimple(42);
   knot.setSeed(42);
   list<int>::const_iterator j = l.begin();
   for (int i = 0; i < 20; i++, j++)
   {
      int m = rand_integer(1, 100);
      int n = r.rand_integer(1, 100);
      ASSERT_EQ(m, n);
      ASSERT_EQ(m, *j);
   }

   // Now, try saving both states, run different sequences and see if we can restart.

   // save states
   pseudorandom s(r);
   pseudorandom t;
   t = saveRandomState();

   // save what would be results running from this state.
   l.clear();
   for (int i = 0; i < 20; i++)
      l.push_back(r.rand_integer(1, 100));

   // scramble random number generators.
   r.sRandSimple(907);
   knot.setSeed(481);
   for (int i = 0; i < 1000; i++)
   {
      r.rand_integer(1, 100);
      rand_integer(1, 100);
   }

   copyRandomState(t);
   j = l.begin();
   for (int i = 0; i < 20; i++, j++)
   {
      int m = s.rand_integer(1, 100);
      int n = t.rand_integer(1, 100);
      int o = rand_integer(1, 100);
      ASSERT_EQ(m, n);
      ASSERT_EQ(m, o);
      ASSERT_EQ(m, *j);
   }
}

TEST(Bfacf3RegressionTest, Bfacf3Run)
{
   clkCigar square;
   threevector<int> v;

   string targetUnkot("0 1 0 -1 1 0 -1 1 -1 -1 2 -1 -1 3 -1 0 3 -1 0 3 0 0 2 0 -1 2 0 -2 2 0 -2 1 0 -3 1 0 -3 0 0 -2 0 0 -1 0 0 0 0 0");
   clkConformationAsList targetConformation;
   targetConformation.readFromText(targetUnkot);
   targetConformation.getVertex(0, v);
   v *= -1;
   targetConformation.translate(v);
   // Set seed to 42 and apply 100 BFACF moves. The resulting conformation should be targetConformation.

   clkConformationBfacf3 knot0(square), knot1(square);
   knot0.setSeed(42);
   for (int i = 0; i < 100; i++)
      knot0.step();
   clkConformationAsList saved(knot0.getComponent(0));
   clkConformationAsList result0(knot0.getComponent(0));
   result0.getVertex(0, v);
   v *= -1;
   result0.translate(v);

   knot1.setSeed(42);
   for (int i = 0; i < 100; i++)
      knot1.step();
   clkConformationAsList result1(knot1.getComponent(0));
   result1.getVertex(0, v);
   v *= -1;
   result1.translate(v);

   ASSERT_TRUE(knot0.getComponent(0) == knot1.getComponent(0)) << "knot0 and knot1 should be same.";
   ASSERT_TRUE(knot0.getComponent(0) == saved) << "knot0 should not have changed when iterating knot1.";
   ASSERT_TRUE(result0 == targetConformation) << "Both knots should be same as target conformation.";
   ASSERT_TRUE(result0 == result1);

   // Can we save the state of the bfacf3 and restart at this same point?
   clkConformationBfacf3 knot2(targetConformation);
   pseudorandom saveRandom = saveRandomState();
   for (int i = 0; i < 100; i++)
      knot2.step();
   clkConformationAsList secondTargetCoformation(knot2.getComponent(0));
   // try restoring random number generator and run bfacf again
   clkConformationBfacf3 knot3(targetConformation);
   copyRandomState(saveRandom);
   for (int i = 0; i < 100; i++)
      knot3.step();
   clkConformationAsList actualConformation(knot3.getComponent(0));
   secondTargetCoformation.getVertex(0, v);
   v *= -1;
   secondTargetCoformation.translate(v);
   actualConformation.getVertex(0, v);
   v *= -1;
   actualConformation.translate(v);
   ASSERT_TRUE(secondTargetCoformation == actualConformation);
}
