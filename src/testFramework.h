#pragma once

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class unit_test{
protected:
	string test_name, pass_message, fail_message;
public:
	virtual bool execute() = 0;
	void print_pass();
	void print_fail();
	unit_test();
};

class test_suite{
private:
	vector<unit_test*> tests;
public:
	void add_test(unit_test* test);
	void run_suite();
	//void add_test(bool (*test)(), string name);
};

