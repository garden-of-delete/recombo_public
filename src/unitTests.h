#pragma once

class fail: public unit_test{
public:
	bool execute();
	fail();
};

class pass: public unit_test{
public:
	bool execute();
	pass();
};

