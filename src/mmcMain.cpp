#include "mmchain.h"
#include "mmc_config.h"
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
	if (argc <= 1){
		print_usage();
		exit(0);
	}

	MmcConfig config;
	config.input_file = argv[1];
	config.output_file = argv[2];

	for (int i = 3; i < argc; i++){
		if (!strcmp(argv[i], "-zmin")){
			config.z_min = atof(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-recomboParas"))
		{
			string paras(argv[i + 1]);
			if (!config.parse_recombo_paras(paras)) {
				cerr << "Error: unrecognized -recomboParas value '" << paras << "'" << endl;
				return 1;
			}
			i++;
		}
		else if (!strcmp(argv[i], "-zmax")){
			config.z_max = atof(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-q")){
			config.q = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-sr")){
			config.swap_ratio = atof(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-s")){
			config.swap_interval = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-n")){
			config.num_samples = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-c")){
			config.steps_between_samples = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-m")){
			config.num_chains = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-w")){
			config.warmup_steps = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-mode")){
			config.sample_mode = *argv[i + 1];
			i++;
		}
		else if (!strcmp(argv[i], "-seed")){
			config.seed = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-minarc")){
			config.min_arc = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-maxarc")){
			config.max_arc = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-targetlength")){
			config.target_recombo_length = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-bfs")){
			config.block_file_size = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "-t")){
			config.time_limit_hours = atoi(argv[i + 1]);
			i++;
		}
		else if (!strcmp(argv[i], "+s")){
			config.suppress_output = true;
		}
		else{
			cerr << "Error: unrecognized option '" << argv[i] << "'" << endl;
			return 1;
		}
	}

	mmchain mmc;
	mmc.initialize(config);
	mmc.run_mmc();

	return 0;
}
