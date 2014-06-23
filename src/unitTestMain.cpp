#include "testFramework.h";
#include "unitTests.h";

using namespace std;

int main(void){
	//create tests first
	unit_test* test1 = new pass();
	unit_test* test2 = new fail();
	//create test framework, (doesn't work? -r.s.)
	test_suite suite;
	//register? tests
	suite.add_test(test1);
	suite.add_test(test2);
	//run tests?
	suite.run_suite();
	//do the delete
	delete test1;
	delete test2;
}