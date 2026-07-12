#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "clkConformationAsList.h"
#include "legacy.h"
#include "writhe.h"

using namespace std;

class WritheTest : public ::testing::Test
{
protected:
    void initFromNewsud(const string& newsud)
    {
        clkConformationAsList conf(newsud);
        n = conf.size();
        verts = new ivector[n];
        threevector<int> v;
        for (int i = 0; i < n; i++)
        {
            conf.getVertex(i, v);
            set_ivector(verts[i], v.getX(), v.getY(), v.getZ());
        }
    }

    void TearDown() override
    {
        delete[] verts;
        verts = NULL;
    }

    ivector* verts = NULL;
    int n = 0;
};

TEST_F(WritheTest, PlanarSquareHasZeroWrithe)
{
    initFromNewsud("enws");
    double ACN = 0;
    double w = writhe(ACN, n, verts, 0);
    EXPECT_DOUBLE_EQ(0.0, w);
    EXPECT_DOUBLE_EQ(0.0, ACN);
}

TEST_F(WritheTest, PlanarRectangleHasZeroWrithe)
{
    initFromNewsud("eennwwss");
    double ACN = 0;
    double w = writhe(ACN, n, verts, 0);
    EXPECT_DOUBLE_EQ(0.0, w);
    EXPECT_DOUBLE_EQ(0.0, ACN);
}

TEST_F(WritheTest, ThreeDimensionalConformationHasNonzeroACN)
{
    // enwusd: a 6-vertex conformation that extends into z
    initFromNewsud("enwusd");
    double ACN = 0;
    double w = writhe(ACN, n, verts, 0);
    EXPECT_DOUBLE_EQ(0.0, w);
    EXPECT_DOUBLE_EQ(0.5, ACN);
}

TEST_F(WritheTest, TenVertexConformationHasNonzeroWrithe)
{
    // euennwwssd: 10-vertex 3D conformation with nonzero writhe
    initFromNewsud("euennwwssd");
    double ACN = 0;
    double w = writhe(ACN, n, verts, 0);
    EXPECT_DOUBLE_EQ(0.25, w);
    EXPECT_NEAR(0.5320471084, ACN, 1e-9);
}

TEST_F(WritheTest, PlanarSquareRadiusOfGyration)
{
    initFromNewsud("enws");
    double rog = radius_of_gyration(n, verts);
    EXPECT_NEAR(1.0 / sqrt(2.0), rog, 1e-9);
}

TEST_F(WritheTest, PlanarRectangleRadiusOfGyration)
{
    initFromNewsud("eennwwss");
    double rog = radius_of_gyration(n, verts);
    EXPECT_NEAR(sqrt(1.5), rog, 1e-9);
}

TEST_F(WritheTest, ThreeDimensionalRadiusOfGyration)
{
    initFromNewsud("euennwwssd");
    double rog = radius_of_gyration(n, verts);
    EXPECT_NEAR(1.2688577540, rog, 1e-9);
}

TEST(WritheOpenTest, ParallelSameDirectionIsZero)
{
    ivector AB[2], CD[2];
    set_ivector(AB[0], 0, 0, 0);
    set_ivector(AB[1], 1, 0, 0);
    set_ivector(CD[0], 0, 1, 0);
    set_ivector(CD[1], 1, 1, 0);
    EXPECT_DOUBLE_EQ(0.0, writhe_open(AB, CD));
}

TEST(WritheOpenTest, AntiparallelIsZero)
{
    ivector AB[2], CD[2];
    set_ivector(AB[0], 0, 0, 0);
    set_ivector(AB[1], 1, 0, 0);
    set_ivector(CD[0], 1, 1, 0);
    set_ivector(CD[1], 0, 1, 0);
    EXPECT_DOUBLE_EQ(0.0, writhe_open(AB, CD));
}

TEST(WritheOpenTest, PerpendicularEdges)
{
    ivector AB[2], CD[2];
    set_ivector(AB[0], 0, 0, 0);
    set_ivector(AB[1], 1, 0, 0);
    set_ivector(CD[0], 0, 0, 1);
    set_ivector(CD[1], 0, 1, 1);
    EXPECT_NEAR(-1.0 / 12.0, writhe_open(AB, CD), 1e-9);
}

TEST(WritheOpenTest, CrossingAboveIsNegative)
{
    ivector AB[2], CD[2];
    set_ivector(AB[0], -1, 0, 0);
    set_ivector(AB[1], 1, 0, 0);
    set_ivector(CD[0], 0, -1, 1);
    set_ivector(CD[1], 0, 1, 1);
    EXPECT_NEAR(-1.0 / 3.0, writhe_open(AB, CD), 1e-9);
}

TEST(WritheOpenTest, CrossingBelowIsPositive)
{
    ivector AB[2], CD[2];
    set_ivector(AB[0], -1, 0, 0);
    set_ivector(AB[1], 1, 0, 0);
    set_ivector(CD[0], 0, -1, -1);
    set_ivector(CD[1], 0, 1, -1);
    EXPECT_NEAR(1.0 / 3.0, writhe_open(AB, CD), 1e-9);
}
