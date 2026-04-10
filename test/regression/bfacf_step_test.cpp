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
