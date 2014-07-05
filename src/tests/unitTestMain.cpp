#include <testFramework.h>

#include "unitTests.h"

using namespace std;

bool identity()
{
	return 1 == 1;
}

bool arithmetic()
{
	return 1 + 1 == 2;
}

bool testAssertions()
{
	ASSERT(true);
	ASSERT(false);
	
	return true;
}

int main(void)
{
	test_suite suite;

	// tests implemented using subclasses of unit_test
	suite.add_test(new pass());
	suite.add_test(new fail());

	// tests implemented using pointers to functions
	suite.add_test(identity, "identity");
	suite.add_test(arithmetic, "arithmetic", "Arithmetic works.", 
		"Arithmetic doesn't work.");
	
	// tests implemented using assertions
	suite.add_test(testAssertions, "assertions", "Everything's cool.", 
		"First it went right, but then it went wrong.");
	
	suite.run_suite();
}

