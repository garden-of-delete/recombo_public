#include <testFramework.h>

#include "argprocessorTests.h"
#include "randomTests.h"

using namespace std;

int main(void)
{
	test_suite suite;

	suite.add_test(testArgprocessor, "argument processor");

	suite.add_test(testRandom, "raw sequence for pseudorandom number generator", "", 
		"Random number generator fails to produce expected sequence.");
	suite.add_test(testRandomInteger, "pseudorandom integers");

	suite.run_suite();
}

