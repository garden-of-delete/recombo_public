#include <gtest/gtest.h>

#include <string>

#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"

using namespace std;

class BfacfStepTest : public ::testing::Test
{
protected:
   void TearDown() override
   {
      delete knot;
   }

   void init(const string& newsud, int seed)
   {
      clkConformationAsList conf(newsud);
      knot = new clkConformationBfacf3(conf);
      knot->setSeed(seed);
   }

   clkConformationBfacf3* knot = NULL;
};

TEST_F(BfacfStepTest, SquareInitializesFromNewsud)
{
   init("enws", 42);
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
}

TEST_F(BfacfStepTest, RectangleInitializesFromNewsud)
{
   init("eennwwss", 42);
   EXPECT_EQ(knot->writeAsNewsud(), "eennwwss");
}

TEST_F(BfacfStepTest, SquareSeed2Step0GrowsByTwo)
{
   // seed 2: first step grows enws -> enwwse (+2)
   init("enws", 2);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
}

TEST_F(BfacfStepTest, SquareSeed2ThreeStepsShrinks)
{
   // seed 2: step 0: enws -> enwwse (+2)
   //         step 1: enwwse -> enwwse (0)
   //         step 2: enwwse -> enwwse (0)
   //         step 3: enwwse -> enws (-2)
   init("enws", 2);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
}

TEST_F(BfacfStepTest, RectangleSeed42FirstStepGrows)
{
   // seed 42: first step grows eennwwss -> euennwwssd (+2, into z dimension)
   init("eennwwss", 42);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "euennwwssd");
}

TEST_F(BfacfStepTest, RectangleSeed2TrajectoryToSquare)
{
   // seed 2: rectangle shrinks back to a square over 8 steps
   //   step 0: eennwwss -> eennwwss (0)
   //   step 1: eennwwss -> sennwwse (0, different start vertex)
   //   step 2: sennwwse -> sennwwse (0)
   //   step 3: sennwwse -> sennwwse (0)
   //   step 4: sennwwse -> enwwse   (-2)
   //   step 5: enwwse   -> enwwse   (0)
   //   step 6: enwwse   -> enwwse   (0)
   //   step 7: enwwse   -> enwwse   (0)
   //   step 8: enwwse   -> enws     (-2, back to square)
   init("eennwwss", 2);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "eennwwss");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "sennwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "sennwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "sennwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
}

TEST_F(BfacfStepTest, SquareSeed1GrowsIntoZDimension)
{
   // seed 1: steps 0-5 are no-ops, step 6 grows into z
   //         enws -> enwusd (+2)
   init("enws", 1);
   for (int i = 0; i < 6; i++)
      knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwusd");
}

TEST_F(BfacfStepTest, NeverShrinksBelowFourVertices)
{
   init("enws", 2);
   for (int i = 0; i < 200; i++)
   {
      knot->step();
      EXPECT_GE((int)knot->writeAsNewsud().length(), 4) << "step " << i;
   }
}

TEST_F(BfacfStepTest, HigherZFavorsGrowth)
{
   // z=0.3 is above DEFAULT_Z (~0.208), so the knot grows more readily
   // At step 2 this diverges from DEFAULT_Z (grows instead of staying)
   init("enws", 2);
   knot->setZ(0.3);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwssen");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "nwssen");
}

TEST_F(BfacfStepTest, HighZGrowsQuickly)
{
   // z=0.4, seed=2: knot grows rapidly compared to default z
   init("enws", 2);
   knot->setZ(0.4);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "endwuwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "endwuwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "ednwuwse");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "denwuwse");
}

TEST_F(BfacfStepTest, HighZSeed42Trajectory)
{
   // z=0.4, seed=42: different seed, same high z
   init("enws", 42);
   knot->setZ(0.4);
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "ednwsu");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "edwu");
   knot->step();
   EXPECT_EQ(knot->writeAsNewsud(), "edwu");
}

TEST_F(BfacfStepTest, StepQWithQ1MatchesStep)
{
   // stepQ with q=1 should produce the same trajectory as step()
   init("enws", 2);
   knot->stepQ(1, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->stepQ(1, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->stepQ(1, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->stepQ(1, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enws");
}

TEST_F(BfacfStepTest, StepQWithQ2SquareSeed2)
{
   // q=2 changes the length-dependent move probabilities
   init("enws", 2);
   knot->init_Q(DEFAULT_Z, 2);
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwse");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enwwssen");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "nwssen");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "nwssen");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "nwssen");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "nwse");
}

TEST_F(BfacfStepTest, StepQWithQ2RectangleSeed42)
{
   init("eennwwss", 42);
   knot->init_Q(DEFAULT_Z, 2);
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "euennwwssd");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "euenwnwssd");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "euennwwssd");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "euennwwssd");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "eunenwwssd");
   knot->stepQ(2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "enenwwss");
}

TEST_F(BfacfStepTest, StepQMultiStepOverload)
{
   // stepQ(long c, int q, double z) runs c steps
   init("enws", 2);
   knot->init_Q(DEFAULT_Z, 2);
   knot->stepQ(7L, 2, DEFAULT_Z);
   EXPECT_EQ(knot->writeAsNewsud(), "nwse");
}

TEST_F(BfacfStepTest, StepQNeverShrinksBelowFourVertices)
{
   init("enws", 2);
   knot->init_Q(DEFAULT_Z, 2);
   for (int i = 0; i < 200; i++)
   {
      knot->stepQ(2, DEFAULT_Z);
      EXPECT_GE((int)knot->writeAsNewsud().length(), 4) << "step " << i;
   }
}
