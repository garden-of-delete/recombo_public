#include "testFramework.h"

#include <runRecombo.h>
#include <recomboArgprocessor.h>

#include <iostream>

#define EPSILON 10e-10
#define DELTA 10e-4

using namespace std;

class basictestclass{
public:
	// argument processors for testRunRecombo
	recomboArgprocessor argsA, argsB, argsH;
	// strings containing initial coordinates in coord format
	std::string trefoil, sixcat;

	basictestclass() 
	{ 
		cout << "basictestclass" << endl;
		
		trefoil = "-1 0 -1\n0 0 -1\n0 1 -1\n0 1 0\n0 1 1\n0 0 1\n0 -1 1\n0 -1 0\n0 -1 -1\n1 -1 -1\n1 0 -1\n1 1 -1\n1 2 -1\n0 2 -1\n0 2 0\n-1 2 0\n-1 1 0\n-1 0 0\n0 0 0\n1 0 0\n1 -1 0\n1 -2 0\n0 -2 0\n-1 -2 0\n-1 -1 0\n-1 -1 -1";
		sixcat = "1 -1 1\n1 -1 0\n1 -1 -1\n0 -1 -1\n0 0 -1\n0 1 -1\n0 2 -1\n-1 2 -1\n-1 2 -2\n-1 2 -3\n-1 1 -3\n-1 0 -3\n-1 0 -2\n-1 0 -1\n-1 0 0\n-1 0 1\n-1 0 2\n0 0 2\n1 0 2\n1 -1 2\n\n2 -2 0\n1 -2 0\n0 -2 0\n0 -1 0\n0 0 0\n1 0 0\n1 0 -1\n1 0 -2\n0 0 -2\n0 1 -2\n-1 1 -2\n-1 1 -1\n-1 1 0\n-1 1 1\n-1 1 2\n-1 1 3\n0 1 3\n0 0 3\n0 -1 3\n0 -1 2\n0 -1 1\n0 0 1\n1 0 1\n2 0 1\n2 -1 1\n2 -1 0";
		vector<string> arglist;

		argsA.splitForArgv("3_1 3_1A -A -s 42 -n 100 --minarc 57 --maxarc 63 -z 0.20961", arglist);
		argsA.process(arglist);

		arglist.clear();
		argsB.splitForArgv("3_1 3_1B -B -s 42 -n 100 --minarc 57 --maxarc 63 -z 0.20961", arglist);
		argsB.process(arglist);

		arglist.clear();
		argsH.splitForArgv("3_1 3_1H -H -s 42 -n 100 --minarc 57 --maxarc 63 -z 0.20961", arglist);
		argsH.process(arglist);
	}

};

bool process(const string& line, ArgumentProcessor& ap)
{
	vector<string> args;
	ap.splitForArgv(line, args);
	return ap.process(args);
}

#define NUMBER_RECOMO_ARGTESTS 2

bool testRecomboArgs()
{
   recomboArgprocessor a[NUMBER_RECOMO_ARGTESTS];
   process("init", a[0]);
   process("init out", a[1]);

   for(int i = 0; i < NUMBER_RECOMO_ARGTESTS; i++)
   {
      ASSERT(a[0].getInitialConformationFileName() == "init");
   }
   
   ASSERT(a[0].getBeforeFileName() == "init_before.b");
   ASSERT(a[0].getAfterFileName() == "init_after.b");
   ASSERT(a[0].getLogFileName() == "init.log");
   ASSERT(a[0].getStatsFileName() == "init_stats.txt");
   
   return true;
}

class testRunRecomboH: public runRecomboH
{
public:
   testRunRecomboH(recomboArgprocessor& args): runRecomboH(args){}
   
   friend bool testRunRecombo();
};

bool testRunRecombo()
{
	basictestclass btc;
	testRunRecomboH testRunKnot(btc.argsH), testRunLink(btc.argsH);
	stringstream knotStream, linkStream;
	knotStream << btc.trefoil;
	linkStream << btc.sixcat;

	pseudorandom r;

	threevector<int> v;
	clkConformationAsList tmp0, tmp1;

	// because the random number generators are seeded with 42, we know what
	// conformations we can expect after iterating
	clkConformationAsList expectedKnot1, expectedKnot2, expectedKnot3, expectedKnot4;
	expectedKnot1.readFromText("126 127 127 126 128 127 126 129 127 127 129 127 127 130 127 127 130 128 128 130 128 128 130 129 128 129 129 128 128 129 127 128 129 127 127 129 127 127 128 127 128 128 127 128 127 128 128 127 129 128 127 129 128 126 130 128 126 130 129 126 130 130 126 130 131 126 129 131 126 129 131 125 129 132 125 129 133 125 130 133 125 130 133 126 130 133 127 130 133 128 130 132 128 129 132 128 128 132 128 128 131 128 127 131 128 127 132 128 126 132 128 126 131 128 126 130 128 126 129 128 127 129 128 128 129 128 129 129 128 129 128 128 129 127 128 129 127 129 129 126 129 129 126 128 128 126 128 128 126 129 128 125 129 128 124 129 128 123 129 127 123 129 126 123 129 126 123 130 126 124 130 126 124 131 127 124 131 127 124 130 127 125 130 127 125 129 127 125 128 126 125 128 126 125 127 126 126 127");
	expectedKnot2.readFromText("127 128 127 127 129 127 127 129 128 127 129 129 126 129 129 126 128 129 126 127 129 126 126 129 127 126 129 127 126 128 127 127 128 128 127 128 129 127 128 129 128 128 130 128 128 130 128 129 130 129 129 130 129 130 129 129 130 129 129 129 129 128 129 128 128 129 128 129 129 128 129 128 128 130 128 127 130 128 127 130 127 126 130 127 126 130 128 125 130 128 125 129 128 125 128 128 126 128 128 127 128 128 127 128 129 127 128 130 128 128 130 128 127 130 128 127 131 128 126 131 128 126 130 128 126 129 128 125 129 127 125 129 127 125 128 126 125 128 126 126 128 126 127 128 126 127 127 127 127 127 127 127 126 127 128 126 126 128 126 126 128 127");
	expectedKnot3.readFromText("127 128 125 127 129 125 127 129 126 126 129 126 126 129 127 126 128 127 126 127 127 126 127 126 126 127 125 125 127 125 125 126 125 125 126 126 126 126 126 126 126 127 127 126 127 127 127 127 128 127 127 128 126 127 128 126 128 129 126 128 129 126 127 129 127 127 130 127 127 130 127 128 130 126 128 130 126 129 130 127 129 130 128 129 130 129 129 130 130 129 130 130 128 130 131 128 130 131 127 129 131 127 128 131 127 128 131 126 128 130 126 129 130 126 129 130 125 128 130 125 127 130 125 127 130 124 126 130 124 126 129 124 126 128 124 126 128 125 126 128 126 127 128 126 127 128 127 127 128 128 126 128 128 126 128 129 126 129 129 126 129 130 126 129 131 126 128 131 126 128 130 126 127 130 126 126 130 126 126 129 126 125 129 126 125 130 127 125 130 127 125 129 128 125 129 128 125 128 128 125 127 128 125 126 128 126 126 128 127 126 128 127 125 128 128 125 128 128 124 127 128 124");
	expectedKnot4.readFromText("128 128 126 127 128 126 126 128 126 126 127 126 125 127 126 124 127 126 124 126 126 124 126 125 125 126 125 126 126 125 126 125 125 126 125 126 126 126 126 127 126 126 128 126 126 128 126 127 128 126 128 128 127 128 129 127 128 129 127 129 128 127 129 128 127 130 128 128 130 128 129 130 129 129 130 130 129 130 131 129 130 131 130 130 131 130 131 130 130 131 130 130 130 130 130 129 131 130 129 131 131 129 130 131 129 130 131 128 130 130 128 129 130 128 129 131 128 128 131 128 128 131 127 128 130 127 128 130 126 127 130 126 126 130 126 125 130 126 125 130 125 124 130 125 124 129 125 123 129 125 123 128 125 123 127 125 122 127 125 122 128 125 122 128 124 122 127 124 122 127 123 122 127 122 122 128 122 123 128 122 123 128 123 124 128 123 124 129 123 125 129 123 125 129 124 125 129 125 126 129 125 127 129 125 127 128 125 127 127 125 127 127 126 127 127 127 127 128 127 127 129 127 127 129 128 127 129 129 126 129 129 126 128 129 126 128 128 126 127 128 127 127 128 127 126 128 127 126 127 127 125 127 127 124 127 128 124 127 128 124 126 128 124 125 128 125 125 128 126 125 128 126 124 129 126 124 129 126 123 129 127 123 130 127 123 130 128 123 130 128 124 130 128 125 130 127 125 130 126 125 130 126 126 129 126 126 129 127 126 128 127 126");

	ASSERT(testRunKnot.readInitialCoords(knotStream));
	// check the z-value. This will be a running theme.
	ASSERT_DOUBLES_EQUAL(0.20961, testRunKnot.conformation->getZ(), EPSILON);

	// mark the state and do some iterations
	clkConformationAsList knot0(testRunKnot.conformation->getComponent(0));
	testRunKnot.markConformation();
	ASSERT_DOUBLES_EQUAL(0.20961, testRunKnot.conformation->getZ(), EPSILON);
	testRunKnot.conformation->step(1000);
	clkConformationAsList knot1(testRunKnot.conformation->getComponent(0));
	//   cout << knot1 << endl;
	ASSERT(expectedKnot1 == knot1);

	// restore the conformation and see that we have knot0
	testRunKnot.restoreConformation();
	ASSERT_DOUBLES_EQUAL(0.20961, testRunKnot.conformation->getZ(), EPSILON);
	ASSERT(knot0 == testRunKnot.conformation->getComponent(0));

	// do another 1000 iterations. Make sure that we don't get the same result as
	// before
	testRunKnot.conformation->step(1000);
	clkConformationAsList knot2(testRunKnot.conformation->getComponent(0));
	//   cout << knot2 << endl;
	ASSERT(expectedKnot2 == knot2);

	// mark the conformation again, iterate and restore
	testRunKnot.markConformation();
	testRunKnot.conformation->step(1000);
	testRunKnot.restoreConformation();
	ASSERT(knot0 != testRunKnot.conformation->getComponent(0));
	//   cout << endl << knot2 << endl << testRunKnot.conformation->getComponent(0) << endl;
	tmp0 = knot2;
	tmp0.getVertex(0, v);
	v *= -1;
	tmp0.translate(v);
	tmp1 = testRunKnot.conformation->getComponent(0);
	tmp1.getVertex(0, v);
	v *= -1;
	tmp1.translate(v);
	//   ASSERT(knot2 == testRunKnot.conformation->getComponent(0));
	ASSERT(tmp0 == tmp1);

	// shake things up
	testRunKnot.conformation->step(1000);

	// mark the conformation again, but this time save the state of the random
	// number generator
	r = saveRandomState();
	testRunKnot.markConformation();
	clkConformationAsList knot3(testRunKnot.conformation->getComponent(0));
	//   cout << knot3 << endl;
	ASSERT(expectedKnot3 == knot3);
	testRunKnot.conformation->step(1000);
	clkConformationAsList knot4(testRunKnot.conformation->getComponent(0));
	//   cout << endl << knot4 << endl;
	ASSERT(expectedKnot4 == knot4);

	// restore the conformation and the random number generator to the state when
	// marked
	testRunKnot.restoreConformation();
	tmp0 = knot3;
	tmp0.getVertex(0, v);
	v *= -1;
	tmp0.translate(v);
	tmp1 = testRunKnot.conformation->getComponent(0);
	tmp1.getVertex(0, v);
	v *= -1;
	tmp1.translate(v);
	ASSERT(tmp0 == tmp1);
	copyRandomState(r);
	// do 1000 iterations and see that we get the same thing
	testRunKnot.conformation->step(1000);
	tmp0 = knot4;
	tmp0.getVertex(0, v);
	v *= -1;
	tmp0.translate(v);
	tmp1 = testRunKnot.conformation->getComponent(0);
	tmp1.getVertex(0, v);
	v *= -1;
	tmp1.translate(v);
	//   ASSERT(tmp0 == tmp1);
	return true;
}

#define NUMBER_RERECOMO_ARGTESTS 10

bool testRerecomboArgs()
{
   rerecomboArgprocessor c, d;
   // it does not make sense to invoke rerecombo without specifying an input file
   ASSERT(!process("", c));
   // it does not make sense to invoke rerecombo without specifying an output file
   ASSERT(!process("init", d));
   
   rerecomboArgprocessor a[NUMBER_RERECOMO_ARGTESTS];
   ASSERT(process("in out", a[0]));
   ASSERT(process("in out --arc 50", a[1]));
   ASSERT(process("in out --minarc 40 --maxarc 60", a[2]));
   ASSERT(process("in out -s 42", a[3]));
   ASSERT(process("in out --minarc 40 -s 21 --maxarc 60", a[4]));
   
   ASSERT(process("in out -2", a[5]));
   ASSERT(process("in out -2 --arc 50", a[6]));
   ASSERT(process("in out --minarc 40 -2 --maxarc 60", a[7]));
   ASSERT(process("in out -s 42 -2", a[8]));
   ASSERT(process("in out --minarc 40 -s 21 --maxarc 60 -2", a[9]));
   
   ASSERT("in" == a[0].getBeforeFileName());
   ASSERT("out" == a[0].getOutputFileName());
   ASSERT("out.b" == a[0].getAfterFileName());
   ASSERT("out.log" == a[0].getLogFileName());
   ASSERT("out_stats.txt" == a[0].getStatsFileName());
   ASSERT_EQUAL(71, a[0].getMinarc());
   ASSERT_EQUAL(77, a[0].getMaxarc());
   ASSERT(a[0].seedWasSetUsingSystemTimer());
   ASSERT(!a[0].twocomonentlinks());

   ASSERT("in" == a[1].getBeforeFileName());
   ASSERT("out" == a[1].getOutputFileName());
   ASSERT("out.b" == a[1].getAfterFileName());
   ASSERT("out.log" == a[1].getLogFileName());
   ASSERT("out_stats.txt" == a[1].getStatsFileName());
   ASSERT_EQUAL(50, a[1].getMinarc());
   ASSERT_EQUAL(50, a[1].getMaxarc());
   ASSERT(a[1].seedWasSetUsingSystemTimer());
   ASSERT(!a[1].twocomonentlinks());

   ASSERT("in" == a[2].getBeforeFileName());
   ASSERT("out" == a[2].getOutputFileName());
   ASSERT("out.b" == a[2].getAfterFileName());
   ASSERT("out.log" == a[2].getLogFileName());
   ASSERT("out_stats.txt" == a[1].getStatsFileName());
   ASSERT_EQUAL(40, a[2].getMinarc());
   ASSERT_EQUAL(60, a[2].getMaxarc());
   ASSERT(a[2].seedWasSetUsingSystemTimer());
   ASSERT(!a[2].twocomonentlinks());

   ASSERT("in" == a[3].getBeforeFileName());
   ASSERT("out" == a[3].getOutputFileName());
   ASSERT("out.b" == a[3].getAfterFileName());
   ASSERT("out.log" == a[3].getLogFileName());
   ASSERT("out_stats.txt" == a[3].getStatsFileName());
   ASSERT_EQUAL(71, a[3].getMinarc());
   ASSERT_EQUAL(77, a[3].getMaxarc());
   ASSERT(!a[3].seedWasSetUsingSystemTimer());
   ASSERT_EQUAL(42, a[3].getSeed());
   ASSERT(!a[3].twocomonentlinks());

   ASSERT("in" == a[4].getBeforeFileName());
   ASSERT("out" == a[4].getOutputFileName());
   ASSERT("out.b" == a[4].getAfterFileName());
   ASSERT("out.log" == a[4].getLogFileName());
   ASSERT("out_stats.txt" == a[4].getStatsFileName());
   ASSERT_EQUAL(40, a[4].getMinarc());
   ASSERT_EQUAL(60, a[4].getMaxarc());
   ASSERT(!a[4].seedWasSetUsingSystemTimer());
   ASSERT_EQUAL(21, a[4].getSeed());
   ASSERT(!a[4].twocomonentlinks());

   ASSERT("in" == a[5].getBeforeFileName());
   ASSERT("out" == a[5].getOutputFileName());
   ASSERT("out.b" == a[5].getAfterFileName());
   ASSERT("out.log" == a[5].getLogFileName());
   ASSERT("out_stats.txt" == a[5].getStatsFileName());
   ASSERT_EQUAL(71, a[5].getMinarc());
   ASSERT_EQUAL(77, a[5].getMaxarc());
   ASSERT(a[5].seedWasSetUsingSystemTimer());
   ASSERT(a[5].twocomonentlinks());

   ASSERT("in" == a[6].getBeforeFileName());
   ASSERT("out" == a[6].getOutputFileName());
   ASSERT("out.b" == a[6].getAfterFileName());
   ASSERT("out.log" == a[6].getLogFileName());
   ASSERT("out_stats.txt" == a[6].getStatsFileName());
   ASSERT_EQUAL(50, a[6].getMinarc());
   ASSERT_EQUAL(50, a[6].getMaxarc());
   ASSERT(a[6].seedWasSetUsingSystemTimer());
   ASSERT(a[6].twocomonentlinks());

   ASSERT("in" == a[7].getBeforeFileName());
   ASSERT("out" == a[7].getOutputFileName());
   ASSERT("out.b" == a[7].getAfterFileName());
   ASSERT("out.log" == a[7].getLogFileName());
   ASSERT("out_stats.txt" == a[7].getStatsFileName());
   ASSERT_EQUAL(40, a[7].getMinarc());
   ASSERT_EQUAL(60, a[7].getMaxarc());
   ASSERT(a[7].seedWasSetUsingSystemTimer());
   ASSERT(a[7].twocomonentlinks());

   ASSERT("in" == a[8].getBeforeFileName());
   ASSERT("out" == a[8].getOutputFileName());
   ASSERT("out.b" == a[8].getAfterFileName());
   ASSERT("out.log" == a[8].getLogFileName());
   ASSERT("out_stats.txt" == a[8].getStatsFileName());
   ASSERT_EQUAL(71, a[8].getMinarc());
   ASSERT_EQUAL(77, a[8].getMaxarc());
   ASSERT(!a[8].seedWasSetUsingSystemTimer());
   ASSERT_EQUAL(42, a[8].getSeed());
   ASSERT(a[8].twocomonentlinks());

   ASSERT("in" == a[9].getBeforeFileName());
   ASSERT("out" == a[9].getOutputFileName());
   ASSERT("out.b" == a[9].getAfterFileName());
   ASSERT("out.log" == a[9].getLogFileName());
   ASSERT("out_stats.txt" == a[9].getStatsFileName());
   ASSERT_EQUAL(40, a[9].getMinarc());
   ASSERT_EQUAL(60, a[9].getMaxarc());
   ASSERT(!a[9].seedWasSetUsingSystemTimer());
   ASSERT_EQUAL(21, a[9].getSeed());
   ASSERT(a[9].twocomonentlinks());
   
   return true;
}
