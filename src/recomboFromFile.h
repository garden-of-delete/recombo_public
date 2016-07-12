#pragma once

#include "mmchain.h"
#include <algorithm>
#include <string>

class recomboFromFile{
private:
	int min_arc, max_arc, n_components, sampling_mode, block_file_mode, current_block_file_number, seed;
	char read_mode;
	char recombo_orientation;
	bool supress_output, info_mode;
	string* infile_name;
	clkConformationAsList initialComp0, initialComp1;
	clkConformationBfacf3* knot;
	pseudorandom siteSelector;
	ifstream* in;
	ofstream* out, *info_file, *sites_file;
	void do_recombo_knots();
	void do_recombo_knots_all();
	void do_recombo_links();
	void do_recombo_links_all();
	bool read_comp_knots(ifstream* in);
	bool read_comp_links(ifstream* in);
	bool inc_filename();
	string get_current_filename();
public:
	recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int n_components, char read_mode, int sampling_mode, int block_file_mode, bool supress_output, bool info_mode, int seed, char orientation);
	void do_recombo();
	~recomboFromFile();
};

