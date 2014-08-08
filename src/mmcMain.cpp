#include "mmchain.h"

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "mmc input_file [operators/options]" << endl;
	cout << "mmc sets up a composite markov chain sampling process using the provided initial clk file." << endl;
	cout << "mmc operators:" << endl;
	cout << "-zmin\tlowest z-value" << endl;
	cout << "-zmax\thighest z-value" << endl;
	cout << "-q\tq-value" << endl;
	cout << "-sr\tlower bound for swap ratio" << endl;
	cout << "-s\tnumber of BFACF steps between swaps" << endl;
	cout << "-n\tnumber of samples per chain before completion" << endl;
	cout << "-c\tnumber of BFACF steps between samples" << endl;
	cout << "-nchains\tinitial number of CMC chains" << endl;
	cout << "-w\tnumber of BFACF steps to take in each chain to warmup" << endl;
	cout << "-m\tsample mode: ['a' Analyze Only], ['s' Sample Only], ['b' Analyze and Sample], ['f' Filter Samples]" << endl;
	cout << "['f' mode ONLY] -minarc\tsample filtering criteria" << endl;
	cout << "['f' mode ONLY] -maxarc\tsample filtering criteria" << endl;
	cout << "['f' mode ONLY] -targetlength\tsample filtering criteria" << endl;
	cout << "-seed\tset seed for random number generator" << endl;
}

int main(int argc, char* argv[]){
	char* infile;
	double zmin = 0, zmax = 0, sr=0;
	int q = 0, s = 0, n = 0, c = 0, nchains = 0, w = 0, m = 0, seed = 0, minarc = 0, maxarc = 0, targetlength = 0;

	if (argc <= 1){
		print_usage();
		mmchain mmc;
		mmc.initialize();
		mmc.run_mmc();
	}
	else{
		infile = argv[1];
		//outfile = argv[2];

		for (int i = 2; i < argc; i++){
			if (!strcmp(argv[i], "-zmin")){
				zmin = atof(argv[i + 1]);
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
			else if (!strcmp(argv[i], "-nchains")){
				nchains = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-w")){
				w = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-m")){
				m = atoi(argv[i + 1]);
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
			else{
				cout << "unrecognized operator/option. Terminating program...";
				return 0;
			}
		}
		//need to write checks for invalid operating parameters
		mmchain mmc;
		mmc.initialize(infile, zmin, zmax, q, sr, s, n, c, nchains, w, m, seed, minarc, maxarc, targetlength);
		mmc.run_mmc();
	}

}