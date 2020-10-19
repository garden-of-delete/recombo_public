#include "mmchain.h"
using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "mmc input_file initial_output_file_name [options]" << endl;
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
	cout << "---------------------------------------------------------------------------------------------------------------" << endl;
	cout << "When selecting values for minarc and maxarc, notice that if the range is too small, " << endl;
	cout << "it would take long time to run, because the samples that satisfy that minarc and maxarc are rare."<< endl;
}

int main(int argc, char* argv[]){
	bool supress_output = false;
	char* infile,
		*outfile,
		mode;

	string recombo_paras; //used to store recomboParas

	int sequence_type, recombo_type;
	//0:inverted, 1:direct
    // 0: standard,
    // 1: positive,
    // 2: negative,
    // 3: automatic,
    // 4: standard+positive,
    // 5: standard+negative,
    // 6: standard+automatic

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
			else if (!strcmp(argv[i], "-recomboParas"))
			{
				recombo_paras.append(argv[i + 1]);

                if (recombo_paras.length() == 2) {
                    if (recombo_paras[0] == 'i')
                    {
                        sequence_type = 0;
                        if (recombo_paras[1] == 's')
                            recombo_type = 0;
                        else if (recombo_paras[1] == 'p')
                            recombo_type = 1;
                        else if (recombo_paras[1] == 'n')
                            recombo_type = 2;
                        else if (recombo_paras[1] == 'a')
                            recombo_type = 3;
                    }
                    else if (recombo_paras[0] == 'd')
                    {
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
                    if (recombo_paras[0] == 'i')
                    {
                        sequence_type = 0;
                        if (recombo_paras[2] == 'p')
                            recombo_type = 4;
                        else if (recombo_paras[2] == 'n')
                            recombo_type = 5;
                        else if (recombo_paras[2] == 'a')
                            recombo_type = 6;
                    }
                    else if (recombo_paras[0] == 'd')
                    {
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
				cout << "unrecognized option " << argv[i] << ", Terminating program...\n";
				return 0;
			}
		}
		//need to write checks for invalid operating parameters
		mmchain mmc;
		mmc.initialize(infile, outfile, zmin, zmax, q, sr, s, n, c, w, m, mode, seed, minarc, maxarc, targetlength, bfs, t, supress_output, sequence_type, recombo_type);
		mmc.run_mmc();
	}

}
