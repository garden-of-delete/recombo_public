#pragma once

#include <iostream>
#include <list>
#include <fstream>

#include "pseudorandom.h"
#include "autocorr.h"
#include "clk.h"
#include "clkCigar.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "conformationAsList.h"
#include "genericConformation.h"

#define critical_z 0.2134

using namespace std;

class search_data{
public:
	double center, std_dev, std_dev_tol , z/*,autocorr_data, autocorr_error, length_var*/;
	int n;
	//scale tol with current length?

	search_data(): center(0), std_dev(0), std_dev_tol(50), z(0), n(5000)/*,autocorr_data(0), autocorr_error(0), length_var(0)*/{}
};

class Analyzer{
	search_data min, max, guess;
	double init_lower, init_upper;
	int target, q, w, n_components, c_steps;
	clkConformationBfacf3* knot;
	clkConformationAsList initialComp0;
	clkConformationAsList initialComp1;

public:
	bool add_initial_conformation(istream&);
	bool Analyzer::add_initial_conformation_from_file(char* filename);
	void reset();
	bool z_from_length(double target, double mean_tol);
	bool length_from_z(search_data* in, bool probe);
	bool initialize_search(double mean_tol);
	bool check_overlap();
	void write_log(ostream);
	void read_table(istream);
	//void generate_table(char* ifilename, char* ofilename, int min, int max, int tol);
	bool write_table(ofstream);
	//void read_from_file(ifstream&);

	Analyzer(char* filename, double Init_lower, double Init_upper, int warmup, int c, int init_q);
};
