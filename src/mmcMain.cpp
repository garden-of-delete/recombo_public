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
	vector<string> errors;

	if (!strcmp(argv[1], "-config")) {
		if (argc < 3) {
			cerr << "Error: -config requires a file path" << endl;
			return 1;
		}
		if (argc > 3) {
			cerr << "Error: -config and CLI options are mutually exclusive" << endl;
			return 1;
		}
		config = MmcConfig::from_json(argv[2], errors);
	} else {
		config = MmcConfig::from_cli(argc, argv, errors);
	}

	if (!errors.empty()) {
		for (size_t i = 0; i < errors.size(); i++)
			cerr << "Error: " << errors[i] << endl;
		return 1;
	}

	mmchain mmc;
	mmc.initialize(config);
	mmc.run_mmc();

	return 0;
}
