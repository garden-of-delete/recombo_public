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
	cout << "-ncomp\tnumber of components in substrate topology" << endl;
	cout << "-seed\tset seed for reproducable recombination" << endl;
	cout << "--m\tinput file format. use --m b for binary (default). use --m t for plain text." << endl;
	cout << "--bfm\tblock file mode. use --bfm 0 for single file mode (default)." << endl 
		<< "\tuse --bfm 1 and enter input file without .b extension for block file mode." << endl;
	cout << "+s\tsupress status output. For use with shell scripts." << endl;
	cout << "-recomboType:" << endl;
	cout << "    i    : do incoherent(knot->knot) and standard recombo" << endl;
	cout << "    iv   : do incoherent and positive virtual recombo" << endl;
	cout << "    ivn  : do incoherent and negative virtual recombo" << endl;
	cout << "    isv  : do incoherent and both standard and positive virtual recombo" << endl;
	cout << "    isvn : do incoherent and both standard and negative virtual recombo" << endl;
	cout << "    c    : do coherent(knot->link OR link->knot) and standard recombo" << endl;
	cout << "    cv   : do coherent and positive virtual recombo" << endl;
	cout << "    cvn  : do coherent and negative virtual recombo" << endl;
	cout << "    csv  : do coherent and both standard and positive virtual recombo" << endl;
	cout << "    csvn : do coherent and both standard and negative virtual recombo" << endl;
	cout << "    ONLY USE i or c IF KNOT SUBSTRATE IS A LINK!!!!" << endl;
}

int main(int argc, char* argv[]){
	int min_arc = 0, max_arc = 0, ncomp = 0, sampling_mode = 0, block_file_mode = 0, seed = 0;
	char* infile = NULL, *outfile = NULL;
	char read_mode = 0;
	string recombo_paras; //used to store recomboType
	int incoOrco = -1,           //1: incoherent, 0: coherent
		stdOrvir = 1,           //1: standard, 0: virtual, 2: both
		virtualDir = 1;       //1: positive, 0: negative

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
		else if (!strcmp(argv[i], "-recomboType")){
			recombo_paras.append(argv[i + 1]);
			if (recombo_paras.length() == 1){
				if (recombo_paras[0] == 'i')				// i
					incoOrco = 1;
				else										// c
					incoOrco = 0;
			}
			else if (recombo_paras.length() == 2){
				if (recombo_paras[0] == 'i')				// iv
					incoOrco = 1;
				else										// cv
					incoOrco = 0;
				stdOrvir = 0;
			}
			else if (recombo_paras.length() == 3) {
				if (recombo_paras[0] == 'i') {
					incoOrco = 1;
					if (recombo_paras[1] == 'v') {			// ivn
						stdOrvir = 0;
						virtualDir = 0;
					}
					else									// isv
						stdOrvir = 2;
				}
				else {// == 'c'
					incoOrco = 0;
					if (recombo_paras[1] == 'v'){			// cvn
						stdOrvir = 0;
						virtualDir = 0;
					}
					else
						stdOrvir = 2;						// csv
				}
			}
			else {//==4
				if (recombo_paras[0] == 'i')				// isvn
					incoOrco = 1;
				else
					incoOrco = 0;							// csvn
				stdOrvir = 2;
				virtualDir = 0;
			}
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
	//set default read mode if none specified
	if (read_mode == 0){
		read_mode = 'b';
	}
	recomboFromFile recombo(min_arc, max_arc, infile, outfile, ncomp, read_mode, sampling_mode, block_file_mode, supress_output, info_mode, seed, incoOrco, stdOrvir, virtualDir);
	recombo.do_recombo();
	return 0;
}
