#include "testFramework.h";
using namespace std;

unit_test::unit_test(){};

void unit_test::print_pass(){
	cout << pass_message << endl;
}

void unit_test::print_fail(){
	cout << fail_message << endl;
}

void test_suite::add_test(unit_test* test){
	tests.push_back(test);
}

void test_suite::run_suite(){
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
	cout << passed << '/' << tests.size() << "tests passed";
}