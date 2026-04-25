#include <gtest/gtest.h>
#include "test_data_util.h"

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <deque>
#include <list>
#include <cmath>

#include <threevector.h>
#include <genericConformation.h>
#include <clkCigar.h>
#include <clkConformationAsList.h>
#include <clkConformationBfacf3.h>
#include <bfacfProbabilities.h>
#include <bfacfProbabilitiesFromZ.h>
#include <bfacfProbabilitiesFromZFixed.h>
#include <pseudorandom.h>

using namespace std;

#define NUMBER_OF_TREFOILS 5
#define TREFOIL_LENGTH 100
#define EPSILON 1.0e-10
#define DELTA 1.0e-4

class testClkConformationBfacf3 : public clkConformationBfacf3
{
public:
	testClkConformationBfacf3(const clkConformationAsList knot) : clkConformationBfacf3(knot) {};
	friend class testClkConformationBfacf3;
};

// --- Fixtures ---

class ClkTestFixture : public ::testing::Test
{
protected:
	void SetUp() override;
	void TearDown() override;

	vector<double> expectedWr;
	vector<double> expectedAcn;
	vector<double> expectedRog;

	string manyTrefoilString;
	string manyTrefoilCoordString;

	vector<vector<int> > trefoilCoords;
	vector<vector<threevector<int> > > trefoilCoordsAsThreevectors;

	clk* clkTrefoil0;
	clk* clkTrefoil1;
	clk* clkTrefoil2;
	clk* clkTrefoil3;
	clk* clkTrefoil4;
};

class Bfacf3TestFixture : public ::testing::Test
{
protected:
	void SetUp() override {
		square2 = clkConformationAsList(square);
		square2.translate(0, 0, 1);
	}

	clkCigar square;
	clkConformationAsList square2;
};

class BfacfProbsTestFixture : public ::testing::Test
{
protected:
	void SetUp() override {
		clkConformationAsList square1(square);
		knot = new testClkConformationBfacf3(square1);
	}

	void TearDown() override {
		delete knot;
	}

	clkCigar square;
	testClkConformationBfacf3* knot;
};

void ClkTestFixture::TearDown()
{
   delete clkTrefoil0;
   delete clkTrefoil1;
   delete clkTrefoil2;
   delete clkTrefoil3;
   delete clkTrefoil4;
}

void ClkTestFixture::SetUp()
{
   json11::Json data = load_test_data("trefoil_conformations.json");

   expectedAcn.resize(NUMBER_OF_TREFOILS);
   expectedRog.resize(NUMBER_OF_TREFOILS);
   expectedWr.resize(NUMBER_OF_TREFOILS);
   trefoilCoords.assign(NUMBER_OF_TREFOILS, vector<int>(3 * TREFOIL_LENGTH));
   trefoilCoordsAsThreevectors.assign(NUMBER_OF_TREFOILS, vector<threevector<int> >(TREFOIL_LENGTH));

   for (int t = 0; t < NUMBER_OF_TREFOILS; t++) {
      expectedRog[t] = data["expected_rog"][t].number_value();
      expectedWr[t] = data["expected_writhe"][t].number_value();
      expectedAcn[t] = data["expected_acn"][t].number_value();
   }

   manyTrefoilString = data["many_trefoil_string"].string_value();

   for (int t = 0; t < NUMBER_OF_TREFOILS; t++) {
      const json11::Json::array& coords = data["conformations"][t].array_items();
      for (int i = 0; i < 3 * TREFOIL_LENGTH; i++)
         trefoilCoords[t][i] = coords[i].int_value();
      for (int i = 0; i < TREFOIL_LENGTH; i++)
         trefoilCoordsAsThreevectors[t][i].set(coords[3*i].int_value(), coords[3*i+1].int_value(), coords[3*i+2].int_value());
   }

   clk** clkTrefoils[] = {&clkTrefoil0, &clkTrefoil1, &clkTrefoil2, &clkTrefoil3, &clkTrefoil4};
   for (int t = 0; t < NUMBER_OF_TREFOILS; t++) {
      *clkTrefoils[t] = new clkConformationAsList(trefoilCoords[t].data(), TREFOIL_LENGTH);

      stringstream ss;
      for (int i = 0; i < TREFOIL_LENGTH; i++)
         ss << trefoilCoords[t][3*i] << " " << trefoilCoords[t][3*i+1] << " " << trefoilCoords[t][3*i+2] << endl;
      ss << t << endl;
      manyTrefoilCoordString += ss.str();
   }
}

static bool vectorsEqual(const vector<threevector<int> >& u, const vector<threevector<int> >& v)
{
   if (u.size() != v.size()) return false;

   for (size_t i = 0; i < u.size(); i++)
      if (u[i] != v[i]) return false;

   return true;
}

static bool vectorsEqual(const vector<int>& u, const vector<threevector<int> >& v)
{
   if (u.size() != 3 * v.size()) return false;

   for (size_t i = 0; i < v.size(); i++)
      if (v[i] != threevector<int>(u[3 * i], u[3 * i + 1], u[3 * i + 2])) return false;

   return true;
}

// --- ClkTestFixture suite: conformation data and geometry ---

TEST_F(ClkTestFixture, Data)
{
	ASSERT_EQ(expectedAcn.size(), NUMBER_OF_TREFOILS) << "Wrong length for expectedAcn.";

	ASSERT_EQ(expectedRog.size(), NUMBER_OF_TREFOILS) << "Wrong length for expectedWr.";
	ASSERT_EQ(expectedWr.size(), NUMBER_OF_TREFOILS) << "Wrong length for expectedRog.";

	ASSERT_EQ(trefoilCoords.size(), NUMBER_OF_TREFOILS) << "Wrong number of trefoils.";
	for (vector<vector<int> >::const_iterator j = trefoilCoords.begin();
	   j != trefoilCoords.end(); j++)
	{
		ASSERT_EQ(j->size(), 3u * TREFOIL_LENGTH) << "Wrong length of trefoil.";
	}

	ASSERT_EQ(trefoilCoordsAsThreevectors.size(), NUMBER_OF_TREFOILS) << "Wrong number of trefoils as threevectors.";
	for (int j = 0; j < NUMBER_OF_TREFOILS; j++)
	{
		ASSERT_EQ(trefoilCoordsAsThreevectors[j].size(), (size_t)TREFOIL_LENGTH) << "Wrong length of trefoil as threevectors.";

		ASSERT_TRUE(vectorsEqual(trefoilCoords[j], trefoilCoordsAsThreevectors[j])) << "Trefoil as threevectors do not have expected coords.";
	}
}

static void testCopyScalarToVectorHelper(const vector<int>& u)
{
   int length = u.size() / 3;
   vector<threevector<int> > v(length);
   copyScalarToVector(u.begin(), u.end(), v.begin());

   ASSERT_TRUE(vectorsEqual(u, v)) << "Copy does not match original.";

   deque<threevector<int> > dq(length);
   copyScalarToVector(u.begin(), u.end(), dq.begin());

   ASSERT_TRUE(dq.size()) << "Deque copy does not match original in number of vertices size.";
   ASSERT_TRUE(equal(v.begin(), v.end(), dq.begin())) << "Deque copy does not match original.";

   list<threevector<int> > l(length);
   copyScalarToVector(u.begin(), u.end(), l.begin());

   ASSERT_TRUE(l.size()) << "List copy does not match original in number of vertices size.";
   ASSERT_TRUE(equal(v.begin(), v.end(), l.begin())) << "List copy does not match original.";
}

TEST_F(ClkTestFixture, CopyScalarToVector)
{
   for (size_t i = 0; i < trefoilCoords.size(); i++)
      testCopyScalarToVectorHelper(trefoilCoords[i]);
}

static void testGenericOutputHelper(const vector<threevector<int> >& u)
{
   stringstream strforvector, strfordeque, strforlist;
   deque<threevector<int> > dq(u.size());
   list<threevector<int> > l(u.size());
   copy(u.begin(), u.end(), dq.begin());
   copy(u.begin(), u.end(), l.begin());

   outputVertices(u.begin(), u.end(), strforvector, "[", "]", "; ", "(", ")", ", ");
   outputVertices(dq.begin(), dq.end(), strfordeque, "[", "]", "; ", "(", ")", ", ");
   outputVertices(l.begin(), l.end(), strforlist, "[", "]", "; ", "(", ")", ", ");

   string expected(strforvector.str());
   ASSERT_EQ(expected, strfordeque.str());
   ASSERT_EQ(expected, strforlist.str());
}

TEST_F(ClkTestFixture, GenericOutput)
{
   stringstream str1, str2, str3, str4;
   string expected1("1 2 3"), expected2("(1, 2, 3)"),
           expected3("1 2 3 1 1 3"), expected4("[(1, 2, 3); (1, 1, 3)]");
   threevector<int> v[2];
   v[0].set(1, 2, 3);
   v[1].set(1, 1, 3);

   outputVertex(v[0], str1);
   outputVertex(v[0], str2, "(", ")", ", ");
   ASSERT_EQ(expected1, str1.str());
   ASSERT_EQ(expected2, str2.str());

   outputVertices(&v[0], &v[2], str3);
   outputVertices(&v[0], &v[2], str4, "[", "]", "; ", "(", ")", ", ");
   ASSERT_EQ(expected3, str3.str());
   ASSERT_EQ(expected4, str4.str());

   for (size_t i = 0; i < trefoilCoordsAsThreevectors.size(); i++)
      testGenericOutputHelper(trefoilCoordsAsThreevectors[i]);
}

TEST_F(ClkTestFixture, Rog)
{
   for (size_t i = 0; i < trefoilCoordsAsThreevectors.size(); i++)
   {
      double r2 = computeRog(trefoilCoordsAsThreevectors[i].begin(), trefoilCoordsAsThreevectors[i].end());
      double r = sqrt(r2);
      ASSERT_NEAR(expectedRog[i], r, DELTA);
   }
}

TEST_F(ClkTestFixture, WrAcn)
{
   // First test that findLast works correctly.
   ASSERT_EQ(0, findLast(0, 1));
   for (int i = -10; i < 10; i++)
      ASSERT_EQ(i, findLast(i, i + 1));
   for (int i = -10; i < 10; i++)
      for (int j = i; j < i + 15; j++)
         ASSERT_EQ(j, findLast(i, j + 1));

   for (size_t i = 0; i < trefoilCoordsAsThreevectors.size(); i++)
   {
      double Wr, Acn;
      computeWrAcn(trefoilCoordsAsThreevectors[i].begin(), trefoilCoordsAsThreevectors[i].end(), Wr, Acn);
      ASSERT_NEAR(expectedWr[i], Wr, DELTA);
      ASSERT_NEAR(expectedAcn[i], Acn, DELTA);
   }
}

TEST_F(ClkTestFixture, ReadFromText)
{
   stringstream stream;
   stream << manyTrefoilString;
   clkConformationAsList knot;

   ASSERT_TRUE(knot.readFromText(stream)) << "Could not read trefoil0";
   ASSERT_TRUE(knot == *clkTrefoil0) << "trefoil0 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromText(stream)) << "Could not read trefoil1";
   ASSERT_TRUE(knot == *clkTrefoil1) << "trefoil1 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromText(stream)) << "Could not read trefoil2";
   ASSERT_TRUE(knot == *clkTrefoil2) << "trefoil2 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromText(stream)) << "Could not read trefoil3";
   ASSERT_TRUE(knot == *clkTrefoil3) << "trefoil3 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromText(stream)) << "Could not read trefoil4";
   ASSERT_TRUE(knot == *clkTrefoil4) << "trefoil4 does not match expected conformation.";

   ASSERT_FALSE(knot.readFromText(stream)) << "There should only be five conformations to read.";
}

TEST_F(ClkTestFixture, ReadFromCoords)
{
   stringstream stream;
   stream << manyTrefoilCoordString;
   clkConformationAsList knot;

   ASSERT_TRUE(knot.readFromCoords(stream)) << "Could not read trefoil0";
   ASSERT_TRUE(knot == *clkTrefoil0) << "trefoil0 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromCoords(stream)) << "Could not read trefoil1";
   ASSERT_TRUE(knot == *clkTrefoil1) << "trefoil1 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromCoords(stream)) << "Could not read trefoil2";
   ASSERT_TRUE(knot == *clkTrefoil2) << "trefoil2 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromCoords(stream)) << "Could not read trefoil3";
   ASSERT_TRUE(knot == *clkTrefoil3) << "trefoil3 does not match expected conformation.";

   ASSERT_TRUE(knot.readFromCoords(stream)) << "Could not read trefoil4";
   ASSERT_TRUE(knot == *clkTrefoil4) << "trefoil4 does not match expected conformation.";

   ASSERT_FALSE(knot.readFromText(stream)) << "There should only be five conformations to read.";
}

// --- Bfacf3TestFixture suite: BFACF algorithm construction, configuration, and execution ---

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

TEST_F(Bfacf3TestFixture, RandomReset)
{
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

TEST_F(Bfacf3TestFixture, Bfacf3)
{
   threevector<int> v;

   clkConformationBfacf3 knot(square);
   ASSERT_EQ(1, knot.size());
   list<clkConformationAsList> components;
   knot.getComponents(components);
   ASSERT_EQ(1u, components.size());
   clkConformationAsList& comp0 = *components.begin();
   ASSERT_EQ(square.size(), comp0.size());
   ASSERT_TRUE(comp0 == knot.getComponent(0));
   comp0.getVertex(0, v);
   v *= -1;
   comp0.translate(v);
   ASSERT_TRUE(square == comp0);

   clkConformationBfacf3 twoComponentLink(square, square2);
   ASSERT_EQ(2, twoComponentLink.size());
   list<clkConformationAsList> components2;
   twoComponentLink.getComponents(components2);
   ASSERT_EQ(2u, components2.size());
   list<clkConformationAsList>::iterator i = components2.begin();
   clkConformationAsList& comp1 = *i;
   i++;
   clkConformationAsList& comp2 = *i;
   ASSERT_TRUE(comp1 == twoComponentLink.getComponent(0));
   ASSERT_TRUE(comp2 == twoComponentLink.getComponent(1));
   comp1.getVertex(0, v);
   v *= -1;
   comp1.translate(v);
   comp2.translate(v);
   ASSERT_TRUE(square == comp1);
   ASSERT_TRUE(square2 == comp2);
}

static void assertBfacfProbabilitiesEqual(const bfacfProbabilities& p1, const bfacfProbabilities& p2)
{
   ASSERT_NEAR(p1.p0(), p2.p0(), EPSILON) << "p0";
   ASSERT_NEAR(p1.p03p2(), p2.p03p2(), EPSILON) << "p_03p2";
   ASSERT_NEAR(p1.p2p0(), p2.p2p0(), EPSILON) << "p_2p0";
   ASSERT_NEAR(p1.p2p02p2(), p2.p2p02p2(), EPSILON) << "p_2p02p2";
   ASSERT_NEAR(p1.p4p2(), p2.p4p2(), EPSILON) << "p_4p2";
   ASSERT_NEAR(p1.m23p2(), p2.m23p2(), EPSILON) << "p_m23p2";
   ASSERT_NEAR(p1.m2(), p2.m2(), EPSILON) << "p_minus2";
   ASSERT_NEAR(p1.p2(), p2.p2(), EPSILON) << "p_plus2";
}

static void assertBfacfProbabilitiesConsistent(const bfacfProbabilities& p)
{
   ASSERT_NEAR(4.0 * p.p2(), p.p4p2(), EPSILON);
   ASSERT_NEAR(p.p0() + 3.0 * p.p2(), p.p03p2(), EPSILON);
   ASSERT_NEAR(p.m2() + 3.0 * p.p2(), p.m23p2(), EPSILON);
   ASSERT_NEAR(2.0 * p.p0(), p.p2p0(), EPSILON);
   ASSERT_NEAR(2.0 * p.p0() + 2.0 * p.p2(), p.p2p02p2(), EPSILON);
}

TEST_F(Bfacf3TestFixture, Bfacf3SetZ)
{
   double expectedDefaultZ = 0.20815;

   // First try a knot.
   clkConformationBfacf3 knot(square2);
   ASSERT_EQ(expectedDefaultZ, DEFAULT_Z);
   ASSERT_EQ(DEFAULT_Z, knot.getZ());
   ASSERT_EQ(DEFAULT_Z, knot.getComponent(0).getZ());

   // Check that the probabilities are what we expect.
   bfacfProbabilitiesFromZFixed pFromDefaultZ, p17(0.17), p19(0.19), p21(0.21);
   ASSERT_EQ(DEFAULT_Z, pFromDefaultZ.z());
   ASSERT_NEAR(0.88497198909891, pFromDefaultZ.m2(), EPSILON);
   ASSERT_NEAR(0.03834267030036, pFromDefaultZ.p2(), EPSILON);
   ASSERT_NEAR(0.46165732969964, pFromDefaultZ.p0(), EPSILON);
   assertBfacfProbabilitiesConsistent(pFromDefaultZ);

   ASSERT_NEAR(0.92021717125242, p17.m2(), EPSILON);
   ASSERT_NEAR(0.02659427624919, p17.p2(), EPSILON);
   ASSERT_NEAR(0.47340572375081, p17.p0(), EPSILON);
   assertBfacfProbabilitiesConsistent(p17);

   ASSERT_NEAR(0.90228277542182, p19.m2(), EPSILON);
   ASSERT_NEAR(0.03257240819273, p19.p2(), EPSILON);
   ASSERT_NEAR(0.46742759180727, p19.p0(), EPSILON);
   assertBfacfProbabilitiesConsistent(p19);

   ASSERT_NEAR(0.88315817362890, p21.m2(), EPSILON);
   ASSERT_NEAR(0.03894727545703, p21.p2(), EPSILON);
   ASSERT_NEAR(0.46105272454297, p21.p0(), EPSILON);
   assertBfacfProbabilitiesConsistent(p21);

   assertBfacfProbabilitiesEqual(pFromDefaultZ, pFromDefaultZ);
   assertBfacfProbabilitiesEqual(pFromDefaultZ, knot.getComponent(0).getProbabilities());

   // Try resetting z to 0.17 on the link level.
   knot.setZ(0.17);
   ASSERT_EQ(0.17, knot.getZ());
   ASSERT_EQ(0.17, knot.getComponent(0).getZ());
   assertBfacfProbabilitiesEqual(p17, knot.getComponent(0).getProbabilities());

   // Try setting z to 0.21 on the component level.
   knot.getComponent(0).setZ(0.21);
   ASSERT_EQ(0.21, knot.getComponent(0).getZ());
   assertBfacfProbabilitiesEqual(p21, knot.getComponent(0).getProbabilities());

   // Next try similar tests for a two-component link.
   clkConformationBfacf3 link(square, square2);
   ASSERT_EQ(DEFAULT_Z, link.getZ());
   ASSERT_EQ(DEFAULT_Z, link.getComponent(0).getZ());
   ASSERT_EQ(DEFAULT_Z, link.getComponent(1).getZ());
   assertBfacfProbabilitiesEqual(pFromDefaultZ, link.getComponent(0).getProbabilities());
   assertBfacfProbabilitiesEqual(pFromDefaultZ, link.getComponent(1).getProbabilities());

   // Reset z to 0.17 on the link level.
   link.setZ(0.17);
   ASSERT_EQ(0.17, link.getZ());
   ASSERT_EQ(0.17, link.getComponent(0).getZ());
   ASSERT_EQ(0.17, link.getComponent(1).getZ());
   assertBfacfProbabilitiesEqual(p17, link.getComponent(0).getProbabilities());
   assertBfacfProbabilitiesEqual(p17, link.getComponent(1).getProbabilities());

   // Reset z to 0.19 on first component and 0.21 on second component.
   link.getComponent(0).setZ(0.19);
   link.getComponent(1).setZ(0.21);
   ASSERT_EQ(0.19, link.getComponent(0).getZ());
   ASSERT_EQ(0.21, link.getComponent(1).getZ());
   assertBfacfProbabilitiesEqual(p19, link.getComponent(0).getProbabilities());
   assertBfacfProbabilitiesEqual(p21, link.getComponent(1).getProbabilities());
}

TEST_F(Bfacf3TestFixture, Bfacf3Run)
{
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

// --- BfacfProbsTestFixture suite: precomputed probability tables ---

TEST_F(BfacfProbsTestFixture, Precomputed)
{
	double z = .2000, q = 1;
	knot->init_Q(z, q);
	for (int n = 4; n < MAX_PRECOMPUTE_LENGTH; n++){
		ASSERT_EQ(knot->probMap[n].p_plus2, (pow((n + 2), (q - 1))*(z * z)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z));
		ASSERT_EQ(knot->probMap[n].p_minus2, pow(n-2, (q - 1)) / (pow(n-2, (q - 1)) + 3.0*pow(n, q - 1) * z * z));
		ASSERT_EQ(knot->probMap[n].p_0, .5 * (knot->probMap[n].p_plus2 + knot->probMap[n].p_minus2));
	}
	z = .1812; q = 3;
	knot->init_Q(z, q);
	for (int n = 4; n < MAX_PRECOMPUTE_LENGTH; n++) {
        ASSERT_EQ(knot->probMap[n].p_plus2,
               (pow((n + 2), (q - 1)) * (z * z)) / (pow(n, (q - 1)) + 3.0 * pow((n + 2), q - 1) * z * z));
        ASSERT_EQ(knot->probMap[n].p_minus2, pow(n - 2, (q - 1)) / (pow(n - 2, (q - 1)) + 3.0 * pow(n, q - 1) * z * z));
        ASSERT_EQ(knot->probMap[n].p_0, .5 * (knot->probMap[n].p_plus2 + knot->probMap[n].p_minus2));
    }
}

TEST_F(BfacfProbsTestFixture, HandComputed)
{
   json11::Json data = load_test_data("bfacf_probabilities.json");
   const json11::Json::array& cases = data["hand_computed"].array_items();

   for (size_t i = 0; i < cases.size(); i++) {
      double z = cases[i]["z"].number_value();
      int q = cases[i]["q"].int_value();
      int n = cases[i]["n"].int_value();
      knot->init_Q(z, q);
      ASSERT_NEAR(knot->probMap[n].p_plus2, cases[i]["p_plus2"].number_value(), EPSILON)
         << "case " << i << ": z=" << z << " q=" << q << " n=" << n;
      ASSERT_NEAR(knot->probMap[n].p_minus2, cases[i]["p_minus2"].number_value(), EPSILON)
         << "case " << i << ": z=" << z << " q=" << q << " n=" << n;
      ASSERT_NEAR(knot->probMap[n].p_0, cases[i]["p_0"].number_value(), EPSILON)
         << "case " << i << ": z=" << z << " q=" << q << " n=" << n;
   }
}
