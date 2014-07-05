#pragma once

#include <testFramework.h>

class fail: public unit_test
{
public:
	fail();
	~fail();

	bool execute();
};

class pass: public unit_test
{
public:
	pass();
	~pass();

	bool execute();
};

