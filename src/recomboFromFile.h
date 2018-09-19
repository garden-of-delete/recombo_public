#pragma once

//#include "mmchain.h"
#include <algorithm>
#include <string>

#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "pseudorandom.h"

class recomboFromFile{
private:
	int min_arc, max_arc, n_components, sampling_mode, block_file_mode, current_block_file_number, seed, incoOrco, stdOrvir, virtualDir;
	char read_mode;
	bool supress_output, info_mode;
	string* infile_name;
	clkConformationAsList initialComp0, initialComp1;
	clkConformationAsList recomboSites;
	clkConformationBfacf3* knot;
	pseudorandom siteSelector;
	ifstream* in;
	ofstream* out, *info_file, *sites_file;
	void writeSitesFile(clkConformationBfacf3* clk, int site_choice);
	void do_recombo_knots();
	void do_recombo_knots_all();
	void do_recombo_links();
	void do_recombo_links_all();
	bool read_comp_knots(ifstream* in);
	bool read_comp_links(ifstream* in);
	bool inc_filename();
	string get_current_filename();
public:
	recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int n_components, char read_mode, int sampling_mode, int block_file_mode, bool supress_output, bool info_mode, int seed, int inco_Or_co, int std_Or_vir, int virtual_Dir);
	void do_recombo();
	~recomboFromFile();
};

