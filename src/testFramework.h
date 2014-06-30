#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

#define ASSERT(assertion) if(!(assertion)) {cout << " in file " << __FILE__ << ", on line " << __LINE__ ; return false;}

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

