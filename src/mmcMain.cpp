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
	cout << "-n\tnumber of samples per chain before completion" << endl;
	cout << "-c\tnumber of BFACF steps between samples" << endl;
	cout << "-m\tinitial number of CMC chains" << endl;
	cout << "-w\tnumber of BFACF steps to take in each chain to warmup" << endl;
	cout << "-mode\tsample mode: ['a' Analyze Only], ['s' Sample Only], ['b' Analyze and Sample], ['f' Filter Samples] ['e' Eternal Sampling] ['m' Block-Mean Sampling]" << endl;
	cout << "['f' mode ONLY] -minarc\tsample filtering criteria" << endl;
	cout << "['f' mode ONLY] -maxarc\tsample filtering criteria" << endl;
	cout << "['f' mode ONLY] -targetlength\tsample filtering criteria" << endl;
	cout << "-seed\tset seed for random number generator" << endl;
	cout << "-bfs\tnumber of conformations to save per .b file" << endl;
	cout << "+s\tsupress status output. For use with shell scripts." << endl;
}

int main(int argc, char* argv[]){
	bool supress_output = false;
	char* infile, *outfile, mode = 0;
	double zmin = 0, zmax = 0, sr=0;
	int q = 0, s = 0, n = 0, c = 0, m = 0, seed = 0, minarc = 0, maxarc = 0, targetlength = 0, bfs = 0;
	long int w = 0;

	if (argc <= 1){
		print_usage();
		mmchain mmc;
		mmc.initialize();
		mmc.run_mmc();
	}
	else{
		infile = argv[1];
		outfile = argv[2];

		for (int i = 3; i < argc; i++){
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
		mmc.initialize(infile, outfile, zmin, zmax, q, sr, s, n, c, w, m, mode, seed, minarc, maxarc, targetlength, bfs, supress_output);
		mmc.run_mmc();
	}

}