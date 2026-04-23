#include "bfacf_config.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "  bfacf input_file output_file [options]" << endl;
	cout << "  bfacf -config config.json" << endl;
	cout << endl;
	cout << "Run BFACF steps on a cubic lattice conformation (knot or 2-component link)." << endl;
	cout << "Configuration is specified either via CLI options or a JSON config file," << endl;
	cout << "but not both." << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -config FILE\tload all parameters from JSON config file" << endl;
	cout << "  -z\t\tfugacity parameter (required)" << endl;
	cout << "  -q\t\tq value for stepQ [default: 1]" << endl;
	cout << "  -n\t\tnumber of BFACF steps to perform (required)" << endl;
	cout << "  -k\t\tsave conformation every k steps [default: save only final]" << endl;
	cout << "  -bfs\t\tconformations per output file [default: all in one file]" << endl;
	cout << "  -seed\t\tRNG seed [default: system time]" << endl;
	cout << "  -fmt\t\toutput format: cube, text, newsud [default: cube]" << endl;
	cout << "  +s\t\tsuppress progress output" << endl;
}

void write_conformation(clkConformationBfacf3& knot, int n_components,
	ostream& out, const string& format) {
	if (format == "cube") {
		clkConformationAsList comp0(knot.getComponent(0));
		comp0.writeAsCube(out);
		if (n_components == 2) {
			clkConformationAsList comp1(knot.getComponent(1));
			comp1.writeAsCube(out);
		}
	} else if (format == "text") {
		clkConformationAsList comp0(knot.getComponent(0));
		comp0.writeAsText(out);
		if (n_components == 2) {
			clkConformationAsList comp1(knot.getComponent(1));
			comp1.writeAsText(out);
		}
	} else if (format == "newsud") {
		out << knot.writeAsNewsud(0) << endl;
		if (n_components == 2)
			out << knot.writeAsNewsud(1) << endl;
	}
}

int main(int argc, char* argv[]){
	if (argc <= 1){
		print_usage();
		return 0;
	}

	BfacfConfig config;

	if (!strcmp(argv[1], "-config")) {
		if (argc < 3) {
			cerr << "Error: -config requires a file path" << endl;
			return 1;
		}
		if (argc > 3) {
			cerr << "Error: -config and CLI options are mutually exclusive" << endl;
			return 1;
		}

		vector<string> json_errors;
		config = BfacfConfig::from_json(argv[2], json_errors);
		if (!json_errors.empty()) {
			for (size_t i = 0; i < json_errors.size(); i++)
				cerr << "Error: " << json_errors[i] << endl;
			return 1;
		}
	} else {
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
			else if (!strcmp(argv[i], "-z")){
				config.z = atof(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-q")){
				config.q = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-n")){
				config.steps = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-k")){
				config.save_interval = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-bfs")){
				config.block_file_size = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-seed")){
				config.seed = atoi(argv[i + 1]);
				i++;
			}
			else if (!strcmp(argv[i], "-fmt")){
				config.output_format = argv[i + 1];
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

	// Load initial conformation
	clkConformationAsList comp0, comp1;
	int n_components = 0;

	ifstream infile(config.input_file.c_str());
	if (!infile) {
		cerr << "Error: cannot open input file: " << config.input_file << endl;
		return 1;
	}
	if (!comp0.readFromCoords(infile)) {
		cerr << "Error: failed to read conformation from: " << config.input_file << endl;
		return 1;
	}
	n_components = 1;
	if (comp1.readFromCoords(infile))
		n_components = 2;
	infile.close();

	// Create BFACF3 engine
	clkConformationBfacf3* knot;
	if (n_components == 2)
		knot = new clkConformationBfacf3(comp0, comp1);
	else
		knot = new clkConformationBfacf3(comp0);

	knot->setSeed(seed);
	knot->setZ(config.z);

	// Log config
	cout << "input_file= " << config.input_file << endl;
	cout << "output_file= " << config.output_file << endl;
	cout << "n_components= " << n_components << endl;
	cout << "z= " << config.z << endl;
	cout << "q= " << config.q << endl;
	cout << "steps= " << config.steps << endl;
	cout << "save_interval= " << config.save_interval << endl;
	cout << "seed= " << seed << endl;
	cout << "output_format= " << config.output_format << endl;
	cout << endl;

	// Determine output mode
	bool binary_output = (config.output_format == "cube");
	int save_interval = config.save_interval > 0 ? config.save_interval : config.steps;
	int block_file_size = config.block_file_size;

	// Output file management
	ofstream out;
	int file_number = 0;
	int conformations_in_file = 0;
	int total_saved = 0;

	auto open_next_file = [&]() {
		if (out.is_open())
			out.close();
		file_number++;
		stringstream ss;
		ss << config.output_file << "n" << file_number;
		if (binary_output)
			ss << ".b";
		else
			ss << ".txt";

		ios_base::openmode mode = ios::out;
		if (binary_output)
			mode |= ios::binary;
		out.open(ss.str().c_str(), mode);

		if (!out) {
			cerr << "Error: cannot open output file: " << ss.str() << endl;
			exit(1);
		}
		conformations_in_file = 0;
	};

	open_next_file();

	// Run BFACF steps
	for (int step = 1; step <= config.steps; step++) {
		knot->stepQ(1, config.q, config.z);

		if (step % save_interval == 0) {
			// Check if we need a new block file
			if (block_file_size > 0 && conformations_in_file >= block_file_size)
				open_next_file();

			write_conformation(*knot, n_components, out, config.output_format);
			conformations_in_file++;
			total_saved++;

			if (!config.suppress_output)
				cout << "step " << step << " / " << config.steps
					<< "  length=" << knot->getComponent(0).size()
					<< "  saved=" << total_saved << endl;
		}
	}

	out.close();

	if (!config.suppress_output)
		cout << endl << "Done. " << total_saved << " conformation(s) saved." << endl;

	delete knot;
	return 0;
}
