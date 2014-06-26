#include "recomboFromFile.h"
#include <string.h>

void print_usage(){
	cout << "Usage:" << endl;
	cout << "recomboFromFile input_file output_file [operators/options]" << endl;
	cout << "recomboFromFile takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a file of conformations " << 
		"after performing recombination with the user specified criteria." << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-min\tminimum arclength" << endl;
	cout << "-max\tmaximum arclength" << endl;
	cout << "-ncomp\tnumber of components" << endl;
}

int main(int argc, char* argv[]){
	int min_arc = 0, max_arc = 0, ncomp = 0; 
	char* infile = NULL, *outfile = NULL;
	
	if (argc < 6){
		print_usage();
		return 0;
	}
	
	infile = argv[2];
	outfile = argv[3];

	for (int i=1; i < argc; i++){
		if (!strcmp(argv[i], "-min")){
			min_arc = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-max")){
			max_arc = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "-ncomp")){
			ncomp = atoi(argv[i+1]);
			i++;
		}
		else{
			cout << "unrecognized operator/option. Terminating program...";
			return 0;
		}
	}

	//sanity check
	if (min_arc > max_arc){
		cout << "mininum arclength must be strictly <= maximum arclength. Terminating program...";
		return 0;
	}

	recomboFromFile recombo(min_arc, max_arc, infile, outfile, ncomp);
	recombo.do_recombo();
	return 0;
}