#include "testFramework.h"
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

int main(void)
{
	test_suite suite;

	suite.add_test(new pass());
	suite.add_test(new fail());

	suite.add_test(identity, "identity");
	suite.add_test(arithmetic, "arithmetic", "Arithmetic works.", "Arithmetic doesn't work.");

	suite.run_suite();
}

