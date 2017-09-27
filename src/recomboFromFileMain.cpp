#include "recomboFromFile.h"
#include <string.h>

void print_usage(){
	cout << "Test Usage:" << endl;
	cout << "recomboFromFile input_file output_file [operators/options]" << endl;
	cout << "recomboFromFile takes an input file of cubic lattice conformations " << 
		"written in Rob Schaerin's CUBE format and outputs a file of conformations " << 
		"after performing recombination with the user specified criteria." << endl << endl;
	cout << "recomboFromFile operators:" << endl;
	cout << "-min\tminimum arclength" << endl;
	cout << "-max\tmaximum arclength" << endl;
	cout << "-ncomp\tnumber of components" << endl;
	cout << "-seed\tset seed for reproducable recombination" << endl;
	cout << "--m\tinput file format. use --m b for binary (default). use --m t for plain text." << endl;
	cout << "--bfm\tblock file mode. use --bfm 0 for single file mode (default)." << endl 
		<< "\tuse --bfm 1 and enter input file without .b extension for block file mode." << endl;
	cout << "+s\tsupress status output. For use with shell scripts." << endl;
	cout << "--orientation ['p' parallel],['a' antiparallel],['u' a or p]" << endl;
}

int main(int argc, char* argv[]){
	int min_arc = 0, max_arc = 0, ncomp = 0, sampling_mode = 0, block_file_mode = 0, seed = 0;
	char* infile = NULL, *outfile = NULL;
	char read_mode = 0;
	char recombo_orientation = 0;
	bool supress_output = false, info_mode = false;
	
	if (argc < 6){
		print_usage();
		return 0;
	}
	
	infile = argv[1];
	outfile = argv[2];

	for (int i=3; i < argc; i++){
		if (!strcmp(argv[i], "-min")){
			min_arc = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "--orientation")){
                        recombo_orientation = *(argv[i+1]);
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
		else if (!strcmp(argv[i], "-seed")){
			seed = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "--m")){
			read_mode = *(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "--ss")){
			sampling_mode = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "--bfm")){
			block_file_mode = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "+s")){
			supress_output = true;
			i++;
		}
		else if (!strcmp(argv[i], "+i")){ //obsolete
			info_mode = true;
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
	//set defauly read mode, orientation if none specified
	if (read_mode == 0){
		read_mode = 'b';
	}
	if (recombo_orientation == 0){
		recombo_orientation = 'a';
	}
	recomboFromFile recombo(min_arc, max_arc, infile, outfile, ncomp, read_mode, sampling_mode, block_file_mode, supress_output, info_mode, seed, recombo_orientation);
	recombo.do_recombo();
	return 0;
}
