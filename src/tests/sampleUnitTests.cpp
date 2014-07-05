#include "sampleUnitTests.h"

#include <iostream>

using namespace std;

pass::pass():
	unit_test("pass", "This is the pass message.", "This is the fail message.")
{}

pass::~pass()
{
	cout << "Delete pass." << endl;
}

bool pass::execute()
{
	return true;
}

fail::fail():
	unit_test("fail", "This is the pass message.", "This is the fail message.")
{}

fail::~fail()
{
	cout << "Delete fail." << endl;
}

bool fail::execute()
{
	return false;
}

