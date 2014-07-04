#include "zAnalyzer.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "zAnalyzer input_file output_file [operators/options]" << endl;
	cout << "zAnalyzer takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a table of average length / Z " << 
		"vaue associations with the provided tolerance" << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-min\tshortest length to find" << endl;
	cout << "-max\tlongest length to find" << endl;
	cout << "-tol\ttolerated error in avg length of the returned z-value" << endl;
}

int main(int argc, char* argv[]){
	Analyzer analyzer;
	char* infile = NULL, *outfile = NULL;
	int min = 0, max = 0, tol = 0;

	if  (argc < 6){
		print_usage();
		return 0;
	}

	infile = argv[1];
	outfile = argv[2];

		for (int i=3; i < argc; i++){
		if (!strcmp(argv[i], "-min")){
			min = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-max")){
			max = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-tol")){
			tol = atoi(argv[i+1]);
			i++;
		}
		else{
			cout << "unrecognized operator/option. Terminating program...";
			return 0;
		}
	}
		//sanity check
	if (min > max){
		cout << "mininum arclength must be strictly <= maximum arclength. Terminating program...";
		return 0;
	}

	analyzer.generate_table(infile, outfile, min, max, tol);
	return 0;
}
