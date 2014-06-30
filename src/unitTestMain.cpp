#include "testFramework.h"
#include "unitTests.h"

using namespace std;

int main(void)
{
	test_suite suite;

	suite.add_test(new pass());
	suite.add_test(new fail());

	suite.run_suite();
}

