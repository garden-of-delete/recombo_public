#include "zAnalyzer.h"
#include "zanalyzer_config.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "  zAnalyzer input_file [options]" << endl;
	cout << "  zAnalyzer -config config.json" << endl;
	cout << endl;
	cout << "zAnalyzer takes an input file of cubic lattice conformations" << endl;
	cout << "written in Rob Scharein's CUBE format and outputs a table of" << endl;
	cout << "average length / Z value associations with the provided tolerance." << endl;
	cout << "Configuration is specified either via CLI options or a JSON config" << endl;
	cout << "file, but not both." << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -config FILE\tload all parameters from JSON config file" << endl;
	cout << "  -l\t\ttarget length" << endl;
	cout << "  -tol\t\ttolerated error in avg length of the returned z-value" << endl;
	cout << "  -zmin\t\tlower bound for the target Z value" << endl;
	cout << "  -zmax\t\tupper bound for the target Z value" << endl;
	cout << "  -c\t\tBFACF steps between samples" << endl;
	cout << "  -w\t\twarmup steps" << endl;
	cout << "  -q\t\tinitial q value" << endl;
	cout << "  -s\t\tset seed (default is current system time)" << endl;
}

int main(int argc, char* argv[]){
	if (argc <= 1){
		print_usage();
		return 0;
	}

	ZAnalyzerConfig config;
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
		config = ZAnalyzerConfig::from_json(argv[2], errors);
	} else {
		config = ZAnalyzerConfig::from_cli(argc, argv, errors);
	}

	if (!errors.empty()) {
		for (size_t i = 0; i < errors.size(); i++)
			cerr << "Error: " << errors[i] << endl;
		return 1;
	}

	// Validate
	vector<string> validation = config.validate();
	if (!validation.empty()) {
		for (size_t i = 0; i < validation.size(); i++)
			cerr << "Error: " << validation[i] << endl;
		return 1;
	}

	// Resolve seed
	int seed = config.seed;
	if (!seed)
		seed = time(NULL);
	srand(seed);

	cout << "seed= " << seed << endl;

	Analyzer analyzer(const_cast<char*>(config.input_file.c_str()),
		config.z_min, config.z_max, config.warmup, config.steps_between_samples, config.q, seed);
	analyzer.z_from_length(config.target_length, config.tolerance);

	return 0;
}
