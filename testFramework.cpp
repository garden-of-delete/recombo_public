#include "testFramework.h";
#include <iostream>;
using namespace std;

bool test_runner::run_test(unit_test& test){
	try{
		test.execute();
		cout <<"Passed" < endl;
	}
	catch(){
		cout << "Failed" << endl;
	}
}