#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

#define ASSERT(assertion) if(!(assertion)) {cout << " in file " << __FILE__ << ", on line " << __LINE__ ; return false;}
#define ASSERT_MESSAGE(message, assertion) if(!(assertion)) {cout << " " << message << " in file " << __FILE__ << ", on line " << __LINE__ ; return false;}
#define ASSERT_EQUAL(x, y) ASSERT((x) == (y))
#define ASSERT_EQUAL_MESSAGE(message, x, y) ASSERT_MESSAGE((message), (x) == (y))
#define ABS(X) ((X) < 0 ? -(X) : (X))
#define ASSERT_DOUBLES_EQUAL(x, y, epsilon) ASSERT(ABS((x) - (y)) < (epsilon))
#define ASSERT_DOUBLES_EQUAL_MESSAGE(message, x, y, epsilon) ASSERT_MESSAGE((message), ABS((x) - (y)) < (epsilon))

class unit_test
{
public:
	unit_test(const string& test_name, const string& pass_message,
		const string& fail_message);
	virtual ~unit_test();

	virtual bool execute() = 0;
	
	void print_name() const;
	void print_pass() const;
	void print_fail() const;

protected:
	const string test_name;
	const string pass_message;
	const string fail_message;
};

class test_suite
{
public:
	test_suite();
	~test_suite();

	void add_test(unit_test* test);
	void add_test(bool (*test)(), const string& test_name,
		const string& pass_message, const string& fail_message);
	void add_test(bool (*test)(), const string& test_name);
	
	void run_suite() const;

private:
	vector<unit_test*> tests;
};

