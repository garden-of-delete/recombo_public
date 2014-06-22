#include "testFramework.h"
#include "unitTests.h"

using namespace std;

pass::pass(){
	pass_message = "this is the pass message";
}

fail::fail(){
	fail_message = "this is the fail message";
}

bool pass::execute(){
	return true;
}

bool fail::execute(){
	return false;
}

