#include "zAnalyzer.h"

void print_usage(){
	cout << "Usage:" << endl;
	cout << "zAnalyzer input_file output_file [operators/options]" << endl;
	cout << "zAnalyzer takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a table of average length / Z " << 
		"vaue associations with the provided tolerance" << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-t\ttarget length" << endl;
	cout << "-e\terror tolerance in length for result" << endl;
}

int main(int argc, char* argv[]){
	Analyzer analyzer();
	char* infile = NULL;

	if  (argc < 5){
		print_usage();
		return 0;
	}

	infile = argv[1];

}