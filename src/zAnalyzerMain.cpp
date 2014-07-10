#include "zAnalyzer.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "zAnalyzer input_file [operators/options]" << endl;
	cout << "zAnalyzer takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a table of average length / Z " << 
		"vaue associations with the provided tolerance" << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-l\ttarget length" << endl;
	cout << "-tol\ttolerated error in avg length of the returned z-value" << endl;
	cout << "-zmin\tbest guess for Z value that is a lower bound for the target Z value" << endl;
	cout << "-zmax\tbest guess for Z valye that is an upper bound for the target Z value" << endl;
	cout << "-c\tbest guess for the number of BFACF steps needed between samples for independence" << endl;
	cout << "-w\tbest guess for an appropriate warmup" << endl;
	cout << "-q\t best guess for an initial q value" << endl;
}

int main(int argc, char* argv[]){
	char* infile = NULL, *outfile = NULL;
	int target_length = 0, tol = 0, warmup = 0, c = 0, q = 0;
	double z_min = 0, z_max = 0;

	if  (argc < 9){
		print_usage();
		return 0;
	}

	infile = argv[1];

		for (int i=3; i < argc; i++){
			if (!strcmp(argv[i], "-l")){
			target_length = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-tol")){
			tol = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-zmax")){
			z_max = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-zmin")){
			z_min = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-q")){
			q = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-c")){
			c = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-w")){
			warmup = atoi(argv[i + 1]);
			i++;
		}
		else{
			cout << "unrecognized operator/option. Terminating program...";
			return 0;
		}
	}
		//sanity check
	if (z_min > z_max){
		cout << "zmin must be strictly greater than zmax. terminating program..." << endl;
		return 0;
	}
	Analyzer analyzer(infile, z_min, z_max, warmup, c, q);
	analyzer.z_from_length(target_length, tol);
	
	return 0;
}
