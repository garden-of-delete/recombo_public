#include <gtest/gtest.h>

#include "clkConformationAsList.h"
#include "clkCigar.h"
#include "threevector.h"

static void expectVertex(const clk& k, int index, int x, int y, int z)
{
   threevector<int> v;
   k.getVertex(index, v);
   EXPECT_EQ(v.getX(), x) << "vertex " << index << " x";
   EXPECT_EQ(v.getY(), y) << "vertex " << index << " y";
   EXPECT_EQ(v.getZ(), z) << "vertex " << index << " z";
}

TEST(NewsudTest, WriteLengthEqualsVertexCount)
{
   clkCigar square;
   clkConformationAsList asList(square);
   EXPECT_EQ(asList.writeAsNewsud().length(), asList.size());

   clkCigar cigar6(6);
   clkConformationAsList asList6(cigar6);
   EXPECT_EQ(asList6.writeAsNewsud().length(), asList6.size());

   clkCigar cigar10(10);
   clkConformationAsList asList10(cigar10);
   EXPECT_EQ(asList10.writeAsNewsud().length(), asList10.size());
}

TEST(NewsudTest, ReadLengthEqualsStringLength)
{
   clkConformationAsList k;
   k.readFromNewsud("enws");
   EXPECT_EQ(k.size(), 4);

   k.readFromNewsud("eennwwss");
   EXPECT_EQ(k.size(), 8);
}

TEST(NewsudTest, WriteSquare)
{
   clkCigar square;
   clkConformationAsList asList(square);
   EXPECT_EQ(asList.writeAsNewsud(), "enws");
}

TEST(NewsudTest, WriteCigar6)
{
   clkCigar cigar(6);
   clkConformationAsList asList(cigar);
   EXPECT_EQ(asList.writeAsNewsud(), "eenwws");
}

TEST(NewsudTest, ReadSquare)
{
   clkConformationAsList k;
   EXPECT_TRUE(k.readFromNewsud("enws"));
   ASSERT_EQ(k.size(), 4);
   expectVertex(k, 0, 0, 0, 0);
   expectVertex(k, 1, 1, 0, 0);
   expectVertex(k, 2, 1, 1, 0);
   expectVertex(k, 3, 0, 1, 0);
}

TEST(NewsudTest, ReadWithStartingCoords)
{
   clkConformationAsList k;
   EXPECT_TRUE(k.readFromNewsud("enws", 5, 10, 3));
   ASSERT_EQ(k.size(), 4);
   expectVertex(k, 0, 5, 10, 3);
   expectVertex(k, 1, 6, 10, 3);
   expectVertex(k, 2, 6, 11, 3);
   expectVertex(k, 3, 5, 11, 3);
}

TEST(NewsudTest, RoundTrip)
{
   clkConformationAsList k;
   k.readFromNewsud("eennwwss");
   std::string result = k.writeAsNewsud();
   EXPECT_EQ(result, "eennwwss");
}

TEST(NewsudTest, RoundTripWithOffset)
{
   clkConformationAsList k;
   k.readFromNewsud("ensw", 3, 4, 5);
   clkConformationAsList k2(k);
   EXPECT_EQ(k2.writeAsNewsud(), "ensw");
}

TEST(NewsudTest, ReadWithUpDown)
{
   clkConformationAsList k;
   EXPECT_TRUE(k.readFromNewsud("eudn"));
   ASSERT_EQ(k.size(), 4);
   expectVertex(k, 0, 0, 0, 0);
   expectVertex(k, 1, 1, 0, 0);
   expectVertex(k, 2, 1, 0, 1);
   expectVertex(k, 3, 1, 0, 0);
}

TEST(NewsudTest, InvalidCharacter)
{
   clkConformationAsList k;
   EXPECT_FALSE(k.readFromNewsud("enxw"));
}

TEST(NewsudTest, EmptyString)
{
   clkConformationAsList k;
   EXPECT_TRUE(k.readFromNewsud(""));
   EXPECT_EQ(k.size(), 1);
}
