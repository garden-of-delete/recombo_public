#include "mmchain.h"
using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "mmc input_file initial_output_file_name [operators/options]" << endl;
	cout << "mmc sets up a composite markov chain sampling process using the provided initial clk file." << endl;
	cout << "mmc operators:" << endl;
	cout << "-zmin\tlowest z-value" << endl;
	cout << "-zmax\thighest z-value" << endl;
	cout << "-q\tq-value" << endl;
	cout << "-sr\tlower bound for swap ratio" << endl;
	cout << "-s\tnumber of BFACF steps between swaps" << endl;
	cout << "[HALTING CRITERIA] -n\tnumber of samples per chain before completion" << endl;
	cout << "[HALTING CRITERIA] -t\tInteger number of hours to run before completion" << endl;
	cout << "-c\tnumber of BFACF steps between samples" << endl;
	cout << "-m\tinitial number of CMC chains. If m = 1, then z-max will be the z-value" << endl;
	cout << "-w\tnumber of BFACF steps to take in each chain to warmup" << endl;
	cout << "-mode\tsample mode: ['a' Analyze Only], ['s' Sample Only], ['b' Analyze and Sample], ['m' Block-Mean Sampling]" << endl;
	cout << "[OPTIONAL] -minarc\tsample filtering criteria" << endl;
	cout << "[OPTIONAL] -maxarc\tsample filtering criteria" << endl;
	cout << "[PREREQ: minarc > 3 and maxarc > minarc] -targetlength\tsample filtering criteria" << endl;
	cout << "-seed\tset seed for random number generator" << endl;
	cout << "-bfs\tnumber of conformations to save per binary block file" << endl;
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
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	cout << "When selecting values for minarc and maxarc, notice that if the range is too small, " << endl;
	cout << "it would take long time to run, because the samples that satisfy that minarc and maxarc are rare."<< endl;
}

int main(int argc, char* argv[]){
	bool supress_output = false;
	char* infile,
		*outfile,
		mode;

	string recombo_paras; //used to store recomboType

	int incoOrco = -1,           //1: incoherent, 0: coherent
		stdOrvir = 1,           //1: standard, 0: virtual, 2: both
		virtualDir = 1;       //1: positive, 0: negative

	double zmin = 0, zmax = 0, sr=0;

	int q = 0, s = 0, n = 0, c = 0, m = 0, seed = 0, minarc = 0, maxarc = 0, targetlength = 0, bfs = 0, t=0;

	long int w = 0;

	if (argc <= 1){
		print_usage();
		exit(0);
	}
	else{
		infile = argv[1];
		outfile = argv[2];

		for (int i = 3; i < argc; i++){
			if (!strcmp(argv[i], "-zmin")){
				zmin = atof(argv[i + 1]);
				i++;
			}
            else if (!strcmp(argv[i], "-recomboType")){
				recombo_paras.append(argv[i + 1]);
				if (recombo_paras.length() == 1){
					if (recombo_paras[0] == 'i')				// i
						incoOrco = 1;
					else 										// c
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
			else if (!strcmp(argv[i], "-zmax")){
				zmax = atof(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-q")){
				q = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-sr")){
				sr = atof(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-s")){
				s = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-n")){
				n = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-c")){
				c = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-m")){
				m = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-w")){
				w = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-mode")){
				mode = *argv[i + 1];
				i++;
			}
			else if (!strcmp(argv[i], "-seed")){
				seed = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-minarc")){
				minarc = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-maxarc")){
				maxarc = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-targetlength")){
				targetlength = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-bfs")){
				bfs = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-t")){
				t = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "+s")){
				supress_output = true;
			}
			else{
				cout << "unrecognized operator/option. Terminating program...";
				return 0;
			}
		}
		//need to write checks for invalid operating parameters
		mmchain mmc;
		mmc.initialize(infile, outfile, zmin, zmax, q, sr, s, n, c, w, m, mode, seed, minarc, maxarc, targetlength, bfs, t, supress_output, incoOrco, stdOrvir, virtualDir);
		mmc.run_mmc();
	}

}
