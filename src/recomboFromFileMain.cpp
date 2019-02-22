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
	cout << "-recomboParas:" << endl;
	cout << "    is    : do incoherent(knot->knot) and standard recombo" << endl;                   //00
	cout << "    ip   : do incoherent and positive virtual recombo" << endl;                        //01
	cout << "    in  : do incoherent and negative virtual recombo" << endl;                         //02
	cout << "    ia  : do incoherent and automatic virtual recombo" << endl;                        //03
	cout << "    isp   : " << endl;                                                                 //04
	cout << "    isn  : " << endl;                                                                  //05
	cout << "    isa  : " << endl;                                                                  //06
	cout << "    ds    : do only coherent(knot->link OR link->knot) standard recombo" << endl;      //10
	cout << "    dp  : do only coherent positive virtual recombo" << endl;                          //11
	cout << "    dn  : do only coherent negative virtual recombo" << endl;                          //12
	cout << "    da  : do only coherent automatic virtual recombo" << endl;                         //13
	cout << "    dsp  : do a mix of coherent standard and positive virtual recombo" << endl;        //14
	cout << "    dsn  : do a mix of coherent standard and negative virtual recombo" << endl;        //15
	cout << "    dsa  : do a mix of coherent standard and automatic virtual recombo" << endl;       //16
	cout << "    ONLY USE is or ds IF KNOT SUBSTRATE IS A LINK!!!!" << endl;
}

int main(int argc, char* argv[]){
	int min_arc = 0, max_arc = 0, ncomp = 0, sampling_mode = 0, block_file_mode = 0, seed = 0;
	char* infile = NULL, *outfile = NULL;
	char read_mode = 0;
	string recombo_paras; //used to store recomboType
	int sequence_type, recombo_type;

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
			if (recombo_paras.length() == 2) {
				if (recombo_paras[0] == 'i') {
					sequence_type = 0;
					if (recombo_paras[1] == 's')
						recombo_type = 0;
					else if (recombo_paras[1] == 'p')
						recombo_type = 1;
					else if (recombo_paras[1] == 'n')
						recombo_type = 2;
					else if (recombo_paras[1] == 'a')
						recombo_type = 3;
				} else if (recombo_paras[0] == 'd') {
					sequence_type = 1;
					if (recombo_paras[1] == 's')
						recombo_type = 0;
					else if (recombo_paras[1] == 'p')
						recombo_type = 1;
					else if (recombo_paras[1] == 'n')
						recombo_type = 2;
					else if (recombo_paras[1] == 'a')
						recombo_type = 3;
				}
			}
			else if (recombo_paras.length() == 3) {
				if (recombo_paras[0] == 'i') {
					sequence_type = 0;
					if (recombo_paras[2] == 'p')
						recombo_type = 4;
					else if (recombo_paras[2] == 'n')
						recombo_type = 5;
					else if (recombo_paras[2] == 'a')
						recombo_type = 6;
				} else if (recombo_paras[0] == 'd') {
					sequence_type = 1;
					if (recombo_paras[2] == 'p')
						recombo_type = 4;
					else if (recombo_paras[2] == 'n')
						recombo_type = 5;
					else if (recombo_paras[2] == 'a')
						recombo_type = 6;
				}
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
	recomboFromFile recombo(min_arc, max_arc, infile, outfile, ncomp, read_mode, sampling_mode, block_file_mode, supress_output, info_mode, seed, sequence_type, recombo_type);
	recombo.do_recombo();
	return 0;
}
