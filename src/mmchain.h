/* 
 * File:   mmchain.h
 * Author: rstolz
 * 
 * Created on April 17th 2013
 * Updated on 6/24/14
 */

#pragma once

#include <iostream>
#include <list>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <string.h>
#include <sstream>

#include "pseudorandom.h"
#include "autocorr.h"
#include "clk.h"
#include "clkCigar.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "conformationAsList.h"
#include "genericConformation.h"

using namespace std;

/**
 * Implimentation of the Multiple Markov Process
 *
 */

class interval_data{
public:
	double delta_b;
	double sucessful_swaps, attempted_swaps;
	interval_data(): delta_b(0), sucessful_swaps(0), attempted_swaps(0){}
};

class chain{
public:
	clkConformationBfacf3* member_knot;
	vector<double> data;
	double z;
	int size; //for debugging purposes

	chain() : member_knot(NULL), z(0), size(0) {};
	chain(const chain& other) : member_knot(other.member_knot), z(other.z), size(other.size), data(other.data) {};
};

class mmchain{
public:
	/**
	* Opens and reads an existing config.txt file. Sets the random number generator's seed directly and calls set_mmc(...).
	*/
	void initialize();

	/*
	* Alternate constructor for use with command line argument processor
	*/
	void initialize(char* in, char* outfile_name, double zmin, double zmax, int q, double sr, int s, int n, int c, long int w, int m, char mode, int seed, int bfs, bool Supress_output);

	/*
	* Constructor for filtering samples to ones that meet recombination criteria
	*/
	void initialize(char* in, char* outfile_name, double zmin, double zmax, int q, double sr, int s, int n, int c, long int w, int m, char mode, int seed, int bfs,
		int Min_arc, int Max_arc, int Target_recombo_length, bool Supress_output);

	/**
	* adds an initial conformation from the given istream. practically speaking, will only be called from outside add_initial_conformation_From_file(...)
	* in debugging and other special scenarios. 
	* @param in the name of the istream containing the coordinate data
	* @return returns false if istream does not contain valid conformation data or data is in the incorrect format
	*/
	bool add_initial_conformation(istream& in);

	/**
	* Calls add_initial_conformation with the provided filename to initialize chains from coordinate file
	* @param filename name of the initial coordinate file (plaintext)
	* @return returns false if the filename can not be opened as an istream, else returns true
	*/
	bool add_initial_conformation_from_file(string& filename);

	/**
	* begins the automated warmup, calibration, and sampling phases. Terminates the process when complete. 
	*/
	void run_mmc();

	/**
	* shows the autocoorelation time for the sampled statistic across all chains. Will be called by run_mmc() if sample mode is 'a'nalyze or 'b'oth. Public for debugging purposes. 
	*/
	void display_results();

private:
	string outfile_name;
	bool supress_output;
	char sample_mode;
	double z_m, z_1, target_swap_ratio;
	clkConformationAsList initialComp0, initialComp1;
	ofstream logfile;
	//when in filtering mode, out2 is for the after recombo conformations
	ofstream out, out2;
	//when in block mean sampling mode, info_file is used to store state information when sampling
	stringstream info_file_stream;
	ofstream info_file;
	//conformationAsList toPrint;
	int seed, m, n_components, swap_interval, n, c, q, min_arc, max_arc, target_recombo_length, block_file_size, block_file_index, current_block_file_number;
	long int w;

	//variables for filter mode only
	long int sample_attempts;

	/**
	* Randomly selects two adjacent chains and attempts a swap using orlandini's swap criteria.
	* Records the outcome of the swap attempt in the appropriate Interval. 
	*/
	void swap();

	/**
	* In the event that no config.txt file is found, this function creates a new config.txt file with some default values, informs the user, and exits 
	*/
	void create_config_file();

	/**
	* COMBO ULTRASETTER for the parameters needed to begin the calibration process. Called by initialize()
	* @param Z_m the largest z value
	* @param Z_1 the smallest z value
	* @param Q the Q value for use in StepQ(...)
	* @param Target_swap_ratio the swap ratio that will be targeted by calibrate_chains(). The number of chains created by calibrate_chains()
	* will approach infinity as this parameter approaches 1. Default value is .8 per M. Szafron. 
	* @param Swap_interval the number of BFACF steps to take between attempting a swap (aka, calling swap())
	* @param Sample_mode controls how the samples are processed. 'a' records the length of the samples into a vector for autocorr analysis. 
	* useful for testing the integrated autocorrelation data of the chains with a given set of operating parameters. 's' samples the conformations
	* directly into the output file and performs no sample length analysis. 'b' performs both 'a' and 's'.
	* @param n number of samples before process termination. 
	* @param c number of steps between sampling attempts.
	* @param m initial number of chains. This value may be increased by calibrate_chains(). The final m value will approach infinity faster than Z_M - Z_1 approaches infinity.
	* @param w number of warmup steps. Per M. Szafron, w needs to be much larger than c. 1 billion is sufficient for trefoils and more complicated topologies. 10+ billion may be required for the unknot / hopf-link.
	* Warmup should is performed with swapping and must complete before calibration can begin.
	*/
	void set_mmc(double Z_m, double Z_1, int Q, double Target_swap_ratio, int Swap_interval, char Sample_mode, int n, int c, int m, long int w);

	/**
	* old purpose: records operating parameters and filenames to log.txt
	* new purpose: records operating parameters to stdout
	* @param buff string buffer containing time information to be printed to cout
	* @return false if logfile could not be created or opened, else returns true
	*/
	void stamp_log(string buff);

	/**
	* called by insert_new_chain(...)
	* allocates a new chain using initialComp0 and initialComp1, themselves already initialized from the initial coordinate file. 
	* sets the z value for the chain and performs a warmup of w steps. 
	* returns a pointer to the newly initialized and warmed up chain to be pushed into vector<chain> chains. 
	* @param z the z value of the new chain
	* @return a pointer to the inialized and warmed up chain
	*/
	chain* alloc_new_chain(double z);

	/**
	* creates a new chain by calling alloc_new_chain(...). Inserts the new chain at the specified position in vector<chain> chains. 
	* creates a new interval and inserts it into the correct position in vector<interval> intervals. 
	* @param i the interval in which to insert the new chain. For example, i = 0 inserts a new chain between chains[0] and chains[1], cooresponding to intervals[0].
	*/
	void insert_new_chain(int i);

	/**
	* Possibly obsolete function. May be useful in the future. Updates interval spacing information using the z values of each chain. 
	*/
	void update_intervals(); 

	/**
	* called by calibrate_chains(). computes confidence intervals for each swap ratio and checks if it is over or includes the target_swap_ratio.
	* calles insert_new_chain(...) at each neccisary interval to create two new invervals, each with a higher swap ratio in the next calibration cycle. 
	* @return true to indicate no new chains were added, false to indicate at least one new chain was created. 
	*/
	bool test_intervals(); 

	/**
	* called by run_mmc(). It will run a simulated sampling session for w steps after the warmup and evaluate the resulting swap ratios for each interval.
	* the process will increase m if neccisary, adding new chains where needed, and repeat this process until a cycle completes with no new chains added. 
	*/
	void calibrate_chains(); 

	/*
	* These three functions are part of Orlandini's described implimentation. They are tested and working already, 
	* and are preserved here in case they are needed in the future. 
	*/
	/*
	void adjust_intervals(); //obsolete
	double find_alpha(); //obsolete
	double check_swap_ratios(double* stable_swap_ratio); //obsolete
	*/

	/**
	* samples according to the specified sample_mode. sample_mode options described in set_mmc(...) documentation.
	*/
	int sample();

	/**
	* finds the size of the conformation stored in each chain. prints the sizes to cout. useful for debugging
	*/
	void update_size();
	
	/**
	* writes conformations to block files of size block_file_size binary clks
	*/
	void write_to_block_file(clkConformationBfacf3* clk);

	/**
	* writes current conformation to block file, performs recombination and writes result to _after blockfile, undoes recombination
	* NOTE: Has a problem with the post-recombination conformations. Needs looking into at some future date. -Stolz
	* @param clk a pointer to the clkconformationBfacf3 object on which to perform recombination
	* @param site_choice an integer value from 0 to number of sites -1 representing the site to perform 
	*
	*/
	
	//deprecated
	//void write_recombination_to_block_file(clkConformationBfacf3* clk, int site_choice);

	//deprecated
	//bool do_recombo_knots(int current_chain);
	//bool do_recombo_links(int current_chain);

	//deprecated
	/*void write_reci_files();*/

	int count_recombo_sites(clkConformationBfacf3* clk);
	int block_mean_sample();

	//private member objects
	vector<interval_data> intervals;
	vector<chain> chains;
};
