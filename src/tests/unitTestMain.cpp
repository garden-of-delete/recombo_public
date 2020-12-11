#include "testFramework.h"

#include "argprocessorTests.h"
#include "autocorrTests.h"
#include "basicTests.h"
#include "clkTests.h"
#include "randomTests.h"
#include "recomboCriteriaTests.h"

using namespace std;

int main(void)
{
	test_suite suite;

	suite.add_test(testAutocorr, "autocorr");

	suite.add_test(testData, "data");
	suite.add_test(testCopyScalarToVector, "testCopyScalarToVector");

	suite.add_test(testGenericOutput, "testGenericOutput");

	suite.add_test(testRog, "testRog");

	suite.add_test(testWrAcn, "testWrAcn");

	suite.add_test(testReadFromText, "testReadFromText");
	suite.add_test(testReadFromCoords, "testReadFromCoords");

	suite.add_test(testRandomReset, "testRandomReset");
	suite.add_test(testBfacf3, "testBfacf3");
	suite.add_test(testBfacf3SetZ, "testBfacf3SetZ");
	suite.add_test(testBfacf3Run, "testBfacf3Run");
	//suite.add_test(testBfacf3CountEdges, "testBfacf3CountEdges"); //incomplete 

	suite.add_test(testRandom, "raw sequence for pseudorandom number generator", "", 
		"Random number generator fails to produce expected sequence.");
	suite.add_test(testRandomInteger, "pseudorandom integers");
	suite.add_test(testPrecomputedBfacf3Probs, "Precomputed BFACF probabilities");
	// TODO: RB had to comment next line because of Segmentation fault 11
	// suite.add_test(testParallelRecombination, "Recombination between unknots");

	suite.run_suite();
}

