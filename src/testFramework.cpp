#include "testFramework.h"

using namespace std;

unit_test::unit_test(const string& test_name, const string& pass_message,
		const string& fail_message):
	test_name(test_name),
	pass_message(pass_message), fail_message(fail_message)
{
	cout << "Base constructor." << endl;
};

unit_test::~unit_test()
{
	cout << "Base destructor." << endl;
};

void unit_test::print_pass()
{
	cout << pass_message << endl;
}

void unit_test::print_fail()
{
	cout << fail_message << endl;
}

test_suite::test_suite(){}

test_suite::~test_suite()
{
	cout << "Delete the tests." << endl;
	for (int i = 0; i < tests.size(); i++)
	{
		delete tests[i];
	}
}

void test_suite::add_test(unit_test* test)
{
	tests.push_back(test);
}

void test_suite::run_suite()
{
	int passed = 0;
	for (int i = 0; i < tests.size(); i++){
		if((tests[i])->execute()){
			(tests[i])->print_pass();
			passed++;
		}
		else{
			(tests[i])->print_fail();
		}
	}
	int total = tests.size();
	int failed = total - passed;
	if(passed == total)
	{
		cout << "Passed all tests out of " << total << '.' << endl; 
	}
	else if(passed == 0)
	{
		cout << "Failed all tests out of " << total  << '.' << endl; 
	}
	else {
		cout << "Passed " << passed
			<<" and failed " << failed 
			<< " out of " << total << " tests." << endl;
	}
}

