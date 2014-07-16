#include "zAnalyzer.h"

Analyzer::Analyzer(char* filename, double Init_lower, double Init_upper, int warmup, int C_steps, int Q){
	add_initial_conformation_from_file(filename);
	init_lower = Init_lower;
	init_upper = Init_upper;
	w = warmup;
	c_steps = C_steps;
	q = Q;
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

   // set z 
   knot->setZ(0);

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
	min.center = min.std_dev /*= min.autocorr_data = min.autocorr_error*/ = 0;
	guess.center = guess.std_dev /*= guess.autocorr_data = guess.autocorr_error*/ = 0;
	max.center = max.std_dev /*= max.autocorr_data = max.autocorr_error = 0*/;
	min.n = guess.n = max.n = 10000;
	min.std_dev_tol = guess.std_dev_tol = /*max.autocorr_error =*/ 25;
}

bool Analyzer::z_from_length(double target_in, double mean_tol){
	int timeout = 20;
	target = target_in;
	cout << endl <<"Initializing Search..." << endl;
	initialize_search(mean_tol);
	cout << endl;
	//while the target is outside the guess window
	while ((guess.center - guess.std_dev > target + mean_tol) || (guess.center + guess.std_dev < target - mean_tol)){
		cout << "=========================================================================" << endl;
		cout << endl << "Step("<< 21 - timeout <<"): Current Z-vals: "<< min.z << "(" << min.center <<")  " << guess.z << "(" << guess.center << ")  " << max.z << "(" << max.center << ")" << endl;
		cout << "=========================================================================" << endl;
		if(target < guess.center){
			max = guess;
			guess.z = exp((log(min.z) + log(guess.z))/2);
			length_from_z(&guess, false);
			check_overlap();
		}

		else if (target > guess.center){
			min = guess;
			guess.z = exp((log(max.z) + log(guess.z)) / 2);
			length_from_z(&guess, false);
			check_overlap();
		}
		timeout--;
		if (!timeout){
			cout << endl << "Search timed out..." << endl;
			return false;
		}
	}
	cout << endl << "Z: " << guess.z << " Q: " << q << "Avg Length: " << guess.center << endl;
	return true;
}

bool Analyzer::length_from_z(search_data* in, bool probe){
	vector<double> data;
	double mean = 0, std_dev = 0;
	//info.error_mean = DBL_MAX;
	int timeout; //hopefully not needed after the code is complete
	//clkConformationBfacf3 knot(*knot);
	knot->setZ(in->z);
	//knot.setSeed(rand() % 20000); //seed 42 used for experimental reproducability
	if (probe == true)
		timeout = 1;
	else
		timeout = 20;
	while(timeout && ((std_dev > in->std_dev_tol) || std_dev == 0)){
		timeout--;
		data.clear();
		cout << "========================================================================="<< endl;
		cout << in->z << " --> ??? w/ q: " << q << endl;
		cout << "Taking " << in->n << " samples every " << c_steps << " iterations." << endl;
		for (int i = 0; i < in->n; i++)		//want to modify to work for two component links
		{
			  if (n_components == 1){
				knot->stepQ(c_steps, q, in->z);
				data.push_back(knot->getComponent(0).size());
			  }
			else if (n_components == 2){
				knot->stepQ(c_steps, q, in->z);
				data.push_back(knot->getComponent(0).size() + knot->getComponent(1).size());
			  }
		}
		//compute mean
		for (int i = 0; i < data.size(); i++){
			mean += data[i];
		}
		mean /= data.size();

		//compute std dev
		for (int i = 0; i < data.size(); i++){
			std_dev += (pow((mean - data[i]), 2) / data.size());
		}
		std_dev = sqrt(std_dev);

		/*info = ac.autocorrelation(data, false);
		in->length_var = ac.computeVariance(data);
		cout << info << endl;
		cout << in->z << " --> " << info.mean << endl;*/
		//if (timeout &&  (info.error_mean > in->std_dev_tol))
			in->n *= 2;
		
	}
	in->center = mean;
	in->std_dev = std_dev;
	/*in->autocorr_data = info.autocorr;
	in->autocorr_error = info.error_autocorr;
	in->center = info.mean;
	in->error = info.error_mean;*/
	cout << in->z << " --> " << in->center << " +- " << in->std_dev << endl;
	return timeout != 0;
}

bool Analyzer::initialize_search(double mean_tol){
	//probe
	reset();
	//double init_lower = .1, init_upper = critical_z;
	double step_size = (init_upper - init_lower)/20; //changed to /20 for unknot;
	//scale with 1/(z-z0(^2)) //pulls
	/*[3/31/13 10:41:09 AM] Reuben Brasher: length, 1/(z-z_0), 1/(z-z_0)^2, 1/(z_z_0)^3...
      [3/31/13 10:44:59 AM] Reuben Brasher: u = 1/(z-z_0)*/
	search_data temp;
	min.z = init_lower;
	max.z = init_upper;
	guess.z = exp((log(max.z) + log(min.z)) / 2);
	cout << "Warming up with " << w << " steps" << endl;
	length_from_z(&min, true);
	knot->stepQ(w, q, max.z);
top:
	length_from_z(&max, true);
	if (max.center < target){
		q += 1;
		goto top;
	}
	length_from_z(&guess, true);
	
/*top:
	length_from_z(&max, true);
	if (max.center < target){
		q += 1;
		reset();
		goto top;
	}
	/*
	while ((max.center - max.std_dev) > target){ 
		temp = max;
		max.z -= step_size;
		knot->init_Q(max.z, q);
		length_from_z(&max, true);
		if (!(max.center - max.std_dev > target)){
			
		}
	}
	if (guess.center == 0){
		length_from_z(&guess, true);
	}
	if (min.center == 0){
		length from
	}*/
	check_overlap();
	//min.std_dev_tol = guess.std_dev_tol = max.std_dev_tol = mean_tol;
	return true;
}

bool Analyzer::check_overlap(){
	int timeout = 10;
	search_data temp;
	cout << "Checking overlap..." << endl;
	while (timeout && (guess.center + guess.std_dev) >= (max.center - max.std_dev)){
		timeout--;
		guess.std_dev_tol /= 2;
		length_from_z(&guess, false);
		if ((guess.center + guess.std_dev) < (max.center - max.std_dev))
			goto part2;
		max.std_dev_tol /= 2;
		length_from_z(&max, false);

	part2:
		while (timeout && (min.center + min.std_dev) >= (guess.center - guess.std_dev)){
			timeout--;
			min.std_dev_tol /= 2;
			length_from_z(&min, false);
			if ((min.center + min.std_dev) < (guess.center - guess.std_dev))
				goto end;
			guess.std_dev_tol /= 2;
			length_from_z(&guess, false);

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
}
	/*
	void Analyzer::generate_table(char* ifilename, char* ofilename, int min, int max, int tol){
	ofstream out;
	ofstream log;
	add_initial_conformation_from_file(ifilename);

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

	}*/