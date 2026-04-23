#include "mmchain.h"
#include "mmc_config.h"
#include <string.h>
#include <stdlib.h>
using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "  mmc input_file output_file [options]" << endl;
	cout << "  mmc -config config.json" << endl;
	cout << endl;
	cout << "mmc sets up a composite markov chain sampling process using the provided" << endl;
	cout << "initial clk file. Configuration is specified either via CLI options or a" << endl;
	cout << "JSON config file, but not both." << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -config FILE\tload all parameters from JSON config file" << endl;
	cout << "  -zmin\t\tlowest z-value" << endl;
	cout << "  -zmax\t\thighest z-value" << endl;
	cout << "  -q\t\tq-value [default: 1]" << endl;
	cout << "  -sr\t\tlower bound for swap ratio [default: 0.8]" << endl;
	cout << "  -s\t\tnumber of BFACF steps between swaps [default: 5]" << endl;
	cout << "  -n\t\tnumber of samples per chain before completion [HALTING]" << endl;
	cout << "  -t\t\tinteger number of hours to run before completion [HALTING]" << endl;
	cout << "  -c\t\tnumber of BFACF steps between samples [default: 20000]" << endl;
	cout << "  -m\t\tinitial number of CMC chains [default: 1]" << endl;
	cout << "  -w\t\tnumber of BFACF steps to warmup [default: 100000]" << endl;
	cout << "  -mode\t\tsample mode: a(nalyze), s(ample), b(oth), m(ean), f(ilter), r(ecombo)" << endl;
	cout << "  -seed\t\tset seed for RNG [default: system time]" << endl;
	cout << "  -bfs\t\tconformations per binary block file" << endl;
	cout << "  -minarc\tmin arclength for reconnection site filtering" << endl;
	cout << "  -maxarc\tmax arclength for reconnection site filtering" << endl;
	cout << "  -targetlength\trequired conformation length (requires -minarc > 3)" << endl;
	cout << "  +s\t\tsuppress status output" << endl;
	cout << endl;
	cout << "  -recomboParas CODE   reconnection parameters:" << endl;
	cout << "    Inverted repeat:  is (standard), ip (positive), in (negative), ia (automatic)" << endl;
	cout << "    Direct repeat:    ds (standard), dp (positive), dn (negative), da (automatic)" << endl;
	cout << "    Mixed:            dsp, dsn, dsa, isp, isn, isa" << endl;
	cout << endl;
	cout << "At least one halting criterion (-n or -t) is required." << endl;
}

int main(int argc, char* argv[]){
	if (argc <= 1){
		print_usage();
		return 0;
	}

	MmcConfig config;

	if (!strcmp(argv[1], "-config")) {
		// JSON config mode
		if (argc < 3) {
			cerr << "Error: -config requires a file path" << endl;
			return 1;
		}
		if (argc > 3) {
			cerr << "Error: -config and CLI options are mutually exclusive" << endl;
			return 1;
		}

		vector<string> json_errors;
		config = MmcConfig::from_json(argv[2], json_errors);
		if (!json_errors.empty()) {
			for (size_t i = 0; i < json_errors.size(); i++)
				cerr << "Error: " << json_errors[i] << endl;
			return 1;
		}
	} else {
		// CLI mode
		if (argc <= 2) {
			print_usage();
			return 0;
		}
		config.input_file = argv[1];
		config.output_file = argv[2];

		for (int i = 3; i < argc; i++){
			if (!strcmp(argv[i], "-config")){
				cerr << "Error: -config and CLI options are mutually exclusive" << endl;
				return 1;
			}
			else if (!strcmp(argv[i], "-zmin")){
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
	}

	mmchain mmc;
	mmc.initialize(config);
	mmc.run_mmc();

	return 0;
}
