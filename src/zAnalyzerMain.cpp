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
		config = ZAnalyzerConfig::from_json(argv[2], json_errors);
		if (!json_errors.empty()) {
			for (size_t i = 0; i < json_errors.size(); i++)
				cerr << "Error: " << json_errors[i] << endl;
			return 1;
		}
	} else {
		// CLI mode
		config.input_file = argv[1];

		for (int i = 2; i < argc; i++){
			if (!strcmp(argv[i], "-config")){
				cerr << "Error: -config and CLI options are mutually exclusive" << endl;
				return 1;
			}
			else if (!strcmp(argv[i], "-l")){
				config.target_length = atoi(argv[i+1]);
				i++;
			}
			else if (!strcmp(argv[i], "-tol")){
				config.tolerance = atoi(argv[i+1]);
				i++;
			}
			else if (!strcmp(argv[i], "-zmax")){
				config.z_max = atof(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-zmin")){
				config.z_min = atof(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-q")){
				config.q = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-c")){
				config.steps_between_samples = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-w")){
				config.warmup = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-s")){
				config.seed = atoi(argv[i + 1]);
				i++;
			}
			else{
				cerr << "Error: unrecognized option '" << argv[i] << "'" << endl;
				return 1;
			}
		}
	}

	// Validate
	vector<string> errors = config.validate();
	if (!errors.empty()) {
		for (size_t i = 0; i < errors.size(); i++)
			cerr << "Error: " << errors[i] << endl;
		return 1;
	}

	// Resolve seed
	int seed = config.seed;
	if (!seed)
		seed = time(NULL);
	srand(seed);

	Analyzer analyzer(const_cast<char*>(config.input_file.c_str()),
		config.z_min, config.z_max, config.warmup, config.steps_between_samples, config.q);
	analyzer.z_from_length(config.target_length, config.tolerance);

	return 0;
}
