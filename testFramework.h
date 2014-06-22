

class unit_test{
public:
	virtual bool execute();
};


class test_runner{
public:
	bool run_test(unit_test& test);
};

