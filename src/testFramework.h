#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class unit_test
{
public:
	unit_test(const string& test_name, const string& pass_message,
		const string& fail_message);
	virtual ~unit_test();

	virtual bool execute() = 0;
	void print_pass();
	void print_fail();

private:
	string test_name, pass_message, fail_message;
};

class simple_test: public unit_test
{
public:
	simple_test(bool (*test)(), const string& test_name,
		const string& pass_message, const string& fail_message);

	virtual bool execute();
	void print_pass();
	void print_fail();
};

class test_suite
{
private:
	vector<unit_test*> tests;
public:
	test_suite();
	~test_suite();

	void add_test(unit_test* test);
	void run_suite();
};

