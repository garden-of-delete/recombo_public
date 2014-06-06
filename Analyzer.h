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
/*#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
//#include "asdf.h"*/

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

void Analyzer::ui(){
	printf("This will be the user interface!");

	return;
}

bool Analyzer::add_initial_conformation(istream& is){
	 if (!initialComp0.readFromCoords(is))
      return false;
   //   cout << initialComp0 << endl;
	 n_components = 1;
   if (initialComp1.readFromCoords(is))
      n_components = 2;

   if (n_components == 2)
      knot = new clkConformationBfacf3(initialComp0, initialComp1);
   else if (n_components == 1)
      knot = new clkConformationBfacf3(initialComp0);

   // set z and seed the random number generator
   knot->setZ(0);
   knot->setSeed(42);
   //knot->setSeed(rand() % 20000);


   return true;
}

bool Analyzer::add_initial_conformation_from_file(char* filename){
	ifstream in;
	in.open(filename);
	if (!in){
		cout << "ERROR: UNABLE TO OPEN FILE";
		return false;
	}
	if (add_initial_conformation(in))
		return true;
	else
		return false;
}

void Analyzer::reset(){
	min.center = min.error = min.autocorr_data = min.autocorr_error = 0;
	guess.center = guess.error = guess.autocorr_data = guess.autocorr_error = 0;
	max.center = max.error = max.autocorr_data = max.autocorr_error = 0;
	min.n = guess.n = max.n = 5000;
	min.error_tol = guess.error_tol = max.autocorr_error = 50;
}

bool Analyzer::z_from_length(double target_in, double mean_tol){
	int timeout = 20;
	type = "length";
	target = target_in;
	cout << endl <<"Initializing Search..." << endl;
	initialize_search(target, mean_tol);
	cout << endl;
	//main while loop
	while ((guess.center + guess.error > target + mean_tol) || guess.center - guess.error < target - mean_tol){

		cout << "=========================================================================" << endl;
		cout << endl << "Step("<< 21 - timeout <<"): Current Z-vals: "<< min.z << "(" << min.center <<")  " << guess.z << "(" << guess.center << ")  " << max.z << "(" << max.center << ")" << endl;
		cout << "=========================================================================" << endl;
		if(target < guess.center){
			max = guess;
			guess.z = (min.z + guess.z)/2;
			length_from_z(&guess, false);
			check_overlap();
		}

		else if (target > guess.center){
			min = guess;
			guess.z = (max.z + guess.z)/2;
			length_from_z(&guess, false);
			check_overlap();
		}
		timeout--;
	}

	if (!timeout)
		return false;
	else
		return true;
}

bool Analyzer::length_from_z(search_data* in, bool probe){
	vector<double> data;
	autocorr ac;
	autocorrInfo info;
	info.autocorr = info.error_autocorr = 10.0;
	//info.error_mean = DBL_MAX;
	int timeout; //hopefully not needed after the code is complete
	//clkConformationBfacf3 knot(*knot);
	knot->setZ(in->z);
	//knot.setSeed(rand() % 20000); //seed 42 used for experimental reproducability
	if (probe == true)
		timeout = 1;
	else
		timeout = 3;
	while(timeout &&  (info.error_mean > in->error_tol)){
		timeout--;
		data.clear();
		cout << "========================================================================="<< endl;
		cout << in->z << " --> ???" << endl;
		cout << "Warming up with 1000000 steps" << endl;
		knot->step(1000000);
		cout << "Taking " << in->n << " samples every " << in->c << " iterations." << endl;
		for (int i = 0; i < in->n; i++)		//want to modify to work for two component links
		{
			  if (n_components == 1){
				knot->step(in->c);
				data.push_back(knot->getComponent(0).size());
			  }
			else if (n_components == 2){
				knot->step(in->c);
				data.push_back(knot->getComponent(0).size() + knot->getComponent(1).size());
			  }
		}
		info = ac.autocorrelation(data, false);
		in->length_var = ac.computeVariance(data);
		cout << info << endl;
		cout << in->z << " --> " << info.mean << endl;
		if (timeout &&  (info.error_mean > in->error_tol))
			in->n *= 2;
	}
	in->autocorr_data = info.autocorr;
	in->autocorr_error = info.error_autocorr;
	in->center = info.mean;
	in->error = info.error_mean;

	return timeout != 0;
}

bool Analyzer::writhe_from_z(search_data* in, bool test){

	return true;
}

bool Analyzer::RoG_from_z(search_data* in, bool test){

	return true;
}

bool Analyzer::initialize_search(double target_length, double mean_tol){
	//probe
	reset();
	//double init_lower = .1, init_upper = critical_z;
	double init_lower = .21508, init_upper = critical_z-.00000000000000001;//for unknot
	double step_size = (init_upper - init_lower)/20; //changed to /20 for unknot;
	//scale with 1/(z-z0(^2)) //pulls
	/*[3/31/13 10:41:09 AM] Reuben Brasher: length, 1/(z-z_0), 1/(z-z_0)^2, 1/(z_z_0)^3...
      [3/31/13 10:44:59 AM] Reuben Brasher: u = 1/(z-z_0)*/
	search_data temp;
	min.z = init_lower;
	max.z = init_upper;

	if (type == "length")
	{
		length_from_z(&min, true);
		while ((min.center + min.error) < target_length){
			temp = min;
			min.z += step_size;
			min.c = 5000;
			length_from_z(&min, true);
		}
		max = min;
		min = temp;
		guess.z = (min.z + max.z)/2;
		length_from_z(&guess, true);
	}
	else if (type == "writhe")
	{

	}
	else if (type == "RoG")
	{

	}
	else
	{
		return false;
	}
	check_overlap();
	min.error_tol = guess.error_tol = max.error_tol = mean_tol;
	return true;
}

bool Analyzer::check_overlap(){
	int timeout = 10;
	search_data temp;
	cout << "Checking overlap..." << endl;
		while (timeout && (guess.center + guess.error) >= (max.center - max.error)){
		timeout--;
		if(type == "length"){
			guess.error_tol /= 2;
			length_from_z(&guess, false);
			if ((guess.center + guess.error) < (max.center - max.error))
				goto part2;
			max.error_tol /= 2;
			length_from_z(&max, false);
		}
		else if (type == "writhe"){
			guess.error_tol /= 2;
			writhe_from_z(&guess, false);
			if ((guess.center + guess.error) < (max.center - max.error))
				continue;
			max.error_tol /= 2;
			writhe_from_z(&max, false);
		}
		else if (type == "RoG"){
			guess.error_tol /= 2;
			RoG_from_z(&guess, false);
			if ((guess.center + guess.error) < (max.center - max.error))
				continue;
			max.error_tol /= 2;
			RoG_from_z(&max, false);
		}
		else
			return false;
		}

		part2:
	while (timeout && (min.center + min.error) >= (guess.center - guess.error)){
		timeout--;
		if(type == "length"){
			min.error_tol /= 2;
			length_from_z(&min, false);
			if ((min.center + min.error) < (guess.center - guess.error))
				goto end;
			guess.error_tol /= 2;
			length_from_z(&guess, false);
		}
		else if (type == "writhe"){
			min.error_tol /= 2;
			writhe_from_z(&min, false);
			if ((min.center + min.error) < (guess.center + guess.error))
				continue;
			guess.error_tol /= 2;
			writhe_from_z(&guess, false);
		}
		else if (type == "RoG"){
			min.error_tol /= 2;
			RoG_from_z(&min, false);
			if ((min.center + min.error) < (guess.center - guess.error))
				continue;
			guess.error_tol /= 2;
			RoG_from_z(&guess, false);
		}
		else
			return false;

		end:
		if (max.center < min.center){
			temp = min;
			min = max;
			max = temp;
		}
		if (min.center > guess.center){
			temp = min;
			min = guess;
			guess = temp;
		}
		if (max.center < guess.center){
			temp = guess;
			guess = max;
			max = temp;
		}
		// if the target goes out of range of the search, purturbs the z values until it comes back into range
		while (min.center > target){
			min.z -= .0025;
			length_from_z(&min, false);
		}
		while (max.center < target){
			max.z += .0025;
			length_from_z(&max, true);
		}
	}
	return true;
}

void Analyzer::generate_table (char* ifilename, char* ofilename, int min, int max, int tol){
	ofstream out;
	ofstream log;
	add_initial_conformation_from_file(ifilename);
	log.open(*ofilename+" log");
	out.open(ofilename);
	log << "TESTING THE LOGGING FUNCTION"<<endl;

	out << type << "," << tol << "," << ifilename << "," << ofilename << endl;
	out << "target,z-val,center,std_deviation,length_var,autocorr,std_deviation,iterations,samples" << endl;
	for (int i = min; i <= max; i += 2)
	{
		z_from_length(i, tol);
		out << i <<","<< guess.z <<"," << guess.center << "," << guess.error << "," <<guess.length_var << "," << guess.autocorr_data << "," << guess.autocorr_error <<","<< guess.c << "," << guess.n << endl;
	}
	/*int i = 0; //code for the roobs
	i = 100;
	z_from_length(i, 1);
	out << i <<","<< guess.z <<"," << guess.center << "," << guess.error << "," << guess.autocorr_data << "," << guess.autocorr_error <<","<< guess.c << "," << guess.n << endl;
	i = 120;
	z_from_length(i, 1);
	out << i <<","<< guess.z <<"," << guess.center << "," << guess.error << "," << guess.autocorr_data << "," << guess.autocorr_error <<","<< guess.c << "," << guess.n << endl;
	i = 150;
	z_from_length(i, 1);
	out << i <<","<< guess.z <<"," << guess.center << "," << guess.error << "," << guess.autocorr_data << "," << guess.autocorr_error <<","<< guess.c << "," << guess.n << endl;
	i = 200;
	z_from_length(i, 5);
	out << i <<","<< guess.z <<"," << guess.center << "," << guess.error << "," << guess.autocorr_data << "," << guess.autocorr_error <<","<< guess.c << "," << guess.n << endl;
	*/
	out.close();
	log.close();
}
