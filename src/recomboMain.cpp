#include "recombo_config.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "pseudorandom.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

void print_usage(){
	cout << "Usage:" << endl;
	cout << "  recombo input_file output_file [options]" << endl;
	cout << "  recombo -config config.json" << endl;
	cout << endl;
	cout << "Perform topological reconnection on a cubic lattice conformation" << endl;
	cout << "(knot or 2-component link). Finds reconnection sites matching arc" << endl;
	cout << "length criteria and performs reconnection at a randomly chosen site." << endl;
	cout << "Configuration is specified either via CLI options or a JSON config" << endl;
	cout << "file, but not both." << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -config FILE\tload all parameters from JSON config file" << endl;
	cout << "  -minarc\tminimum arc length (required)" << endl;
	cout << "  -maxarc\tmaximum arc length (required)" << endl;
	cout << "  -targetlength\tonly reconnect if conformation has this length" << endl;
	cout << "  -seed\t\tRNG seed [default: system time]" << endl;
	cout << "  -fmt\t\toutput format: cube, text, newsud [default: cube]" << endl;
	cout << "  +s\t\tsuppress progress output" << endl;
	cout << endl;
	cout << "  -recomboParas CODE   reconnection parameters (required):" << endl;
	cout << "    Inverted repeat:  is (standard), ip (positive), in (negative), ia (automatic)" << endl;
	cout << "    Direct repeat:    ds (standard), dp (positive), dn (negative), da (automatic)" << endl;
	cout << "    Mixed:            dsp, dsn, dsa, isp, isn, isa" << endl;
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

void write_sites(clkConformationBfacf3& knot, int site_choice, int n_components,
	ostream& out) {
	vector<threevector<int> > sitelist = knot.getChosenSite(site_choice);
	clkConformationAsList sites;

	if ((knot.dirIsParal(site_choice) && n_components == 1) ||
		(knot.dirIsParal(site_choice) && n_components == 2)) {
		sites.addVertexBack(sitelist[0]);
		sites.addVertexBack(sitelist[1]);
		sites.addVertexBack(sitelist[3]);
		sites.addVertexBack(sitelist[2]);
	} else {
		for (int j = 0; j < 4; ++j)
			sites.addVertexBack(sitelist[j]);
	}

	sites.writeAsCube(out);
}

int main(int argc, char* argv[]){
	if (argc <= 1){
		print_usage();
		return 0;
	}

	RecomboConfig config;
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
		config = RecomboConfig::from_json(argv[2], errors);
	} else {
		config = RecomboConfig::from_cli(argc, argv, errors);
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

	// Create BFACF3 engine (needed for recombo site detection)
	clkConformationBfacf3* knot;
	if (n_components == 2)
		knot = new clkConformationBfacf3(comp0, comp1);
	else
		knot = new clkConformationBfacf3(comp0);

	knot->setSeed(seed);

	// Check target length constraint
	if (config.target_length > 0) {
		int length = knot->getComponent(0).size();
		if (n_components == 2)
			length += knot->getComponent(1).size();
		if (length != config.target_length) {
			cerr << "Error: conformation length " << length
				<< " does not match target length " << config.target_length << endl;
			delete knot;
			return 1;
		}
	}

	// Log config
	cout << "input_file= " << config.input_file << endl;
	cout << "output_file= " << config.output_file << endl;
	cout << "n_components= " << n_components << endl;
	cout << "min_arc= " << config.min_arc << endl;
	cout << "max_arc= " << config.max_arc << endl;
	cout << "sequence_type= " << config.sequence_type << endl;
	cout << "recombo_type= " << config.recombo_type << endl;
	cout << "seed= " << seed << endl;
	cout << "output_format= " << config.output_format << endl;
	int length = knot->getComponent(0).size();
	if (n_components == 2)
		length += knot->getComponent(1).size();
	cout << "conformation_length= " << length << endl;
	cout << endl;

	// Find reconnection sites
	int total_para = 0, total_anti = 0, para = 0, anti = 0;
	int sites = knot->countRecomboSites(config.min_arc, config.max_arc,
		config.sequence_type, config.recombo_type,
		total_para, total_anti, para, anti);

	cout << "sites_found= " << sites << endl;
	cout << "parallel_sites= " << para << endl;
	cout << "anti_parallel_sites= " << anti << endl;

	if (sites == 0) {
		cout << "No reconnection sites found. No output written." << endl;
		delete knot;
		return 0;
	}

	// Choose a site randomly
	pseudorandom rng(seed);
	int choice = rng.rand_integer(0, sites);
	cout << "chosen_site= " << choice << endl;
	cout << endl;

	// Open output files
	bool binary = (config.output_format == "cube");
	string ext = binary ? ".b" : ".txt";
	ios_base::openmode mode = ios::out;
	if (binary)
		mode |= ios::binary;

	string before_path = config.output_file + "_before" + ext;
	string after_path = config.output_file + "_after" + ext;
	string sites_path = config.output_file + "_sites.b";

	ofstream before_file(before_path.c_str(), mode);
	if (!before_file) {
		cerr << "Error: cannot open output file: " << before_path << endl;
		delete knot;
		return 1;
	}

	ofstream after_file(after_path.c_str(), mode);
	if (!after_file) {
		cerr << "Error: cannot open output file: " << after_path << endl;
		delete knot;
		return 1;
	}

	ofstream sites_file(sites_path.c_str(), ios::out | ios::binary);
	if (!sites_file) {
		cerr << "Error: cannot open output file: " << sites_path << endl;
		delete knot;
		return 1;
	}

	// Write before conformation
	write_conformation(*knot, n_components, before_file, config.output_format);
	before_file.close();

	// Write sites
	write_sites(*knot, choice, n_components, sites_file);
	sites_file.close();

	// Perform reconnection
	if (knot->performRecombination(cout, config.sequence_type, config.recombo_type, n_components, choice)) {
		// Write after conformation
		write_conformation(*knot, n_components, after_file, config.output_format);
		if (!config.suppress_output)
			cout << "Reconnection performed successfully." << endl;
	} else {
		cerr << "Warning: reconnection failed at chosen site." << endl;
	}

	after_file.close();
	delete knot;
	return 0;
}
