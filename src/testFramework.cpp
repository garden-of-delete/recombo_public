#include "testFramework.h"

#include <vector>

using namespace std;

unit_test::unit_test(const string& test_name, const string& pass_message,
		const string& fail_message):
	test_name(test_name),
	pass_message(pass_message), fail_message(fail_message)
{};

unit_test::~unit_test()
{};

void unit_test::print_name() const
{
	cout << test_name;
}

void print_message(const string& message)
{
	if(!message.empty())
	{
		cout << '\t' << message << endl;
	}
}

void unit_test::print_pass() const
{
	print_message(pass_message);
}

void unit_test::print_fail() const
{
	print_message(fail_message);
}

class simple_test: public unit_test
{
public:
	simple_test(bool (*test)(), const string& test_name,
		const string& pass_message, const string& fail_message):
		unit_test(test_name, pass_message, fail_message),
		test(test)
	{}
	
	simple_test(bool (*test)(), const string& test_name):
		unit_test(test_name, "", ""),
		test(test)
	{}

	~simple_test()
	{
		cout << "Delete " << test_name << '.' << endl;
	}

	bool execute()
	{
		return test();
	}
	
private:
	bool (*test)();
};

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

void test_suite::add_test(bool (*test)(), const string& test_name,
	const string& pass_message, const string& fail_message)
{
	add_test(new simple_test(test, test_name, pass_message, fail_message));
}

void test_suite::add_test(bool (*test)(), const string& test_name)
{
	add_test(new simple_test(test, test_name));
}

void test_suite::run_suite() const
{
	int passed = 0;
	for (int i = 0; i < tests.size(); i++)
	{
		cout << "Test ";
		tests[i]->print_name();
		if((tests[i])->execute()) 
		{
			cout << " passed." << endl;
			tests[i]->print_pass();
			passed++;
		}
		else 
		{
			cout << " failed." << endl;
			tests[i]->print_fail();
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
	else 
	{
		cout << "Passed " << passed
			<<" and failed " << failed 
			<< " out of " << total << " tests." << endl;
	}
}

