#pragma once

#include <iostream>
#include <list>
#include <fstream>

#include "random.h"
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
	double center, error, error_tol, z, autocorr_data, autocorr_error, length_var;
	int n, c;
	//scale tol with current length?

	search_data(): center(0), error(0), error_tol(50), z(0), n(5000), c(2500), autocorr_data(0), autocorr_error(0), length_var(0){}
};

class Analyzer{
	search_data min, max, guess;
	double target;
	char* type;
	clkConformationBfacf3* knot;
	int n_components;
	clkConformationAsList initialComp0;
	clkConformationAsList initialComp1;

	//clk& initial_conformation_2;

public:
	void ui();
	bool add_initial_conformation(istream&);
	bool Analyzer::add_initial_conformation_from_file(char* filename);
	void reset();
	bool z_from_length(double, double);
	bool z_from_writhe(double, double, double*);
	bool z_from_RoG(double, double, double*);
	bool length_from_z(search_data*, bool);
	bool writhe_from_z(search_data*, bool);
	bool RoG_from_z(search_data*, bool);
	bool initialize_search(double, double);
	bool check_overlap();
	void write_log(ostream);
	void read_table(istream);
	void generate_table(char*, char*, int, int, int);
	bool write_table(ofstream);
	//void read_from_file(ifstream&);

	Analyzer (): min(), guess(), max(), type("none"), n_components(0){}
	//Analyzer(clk& init):z(0), n(0), c(0), n_conformations(1), initial_conformation(init), initial_conformation_2(init){}
	//Analyzer(clk& init, clk& init2):z(0), n(0), c(0), n_conformations(2), initial_conformation(init), initial_conformation_2(init2){}
};