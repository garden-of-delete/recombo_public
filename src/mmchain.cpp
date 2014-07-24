#include "mmchain.h"
using namespace std;

void mmchain::initialize(){
	//check if config file exists
	ifstream config("config.txt", ios::in);
	if (config.fail()){
		create_config_file();
		cout << "config.txt not found. It has been created. Enter operational parameters and relaunch mmc.";
		exit(1);
	}
	//read data from config file
	string filename, buff;
	buff.resize(10000000);
	double zmin, zmax, target_swap_ratio, warmup;
	int q, swap_interval, n_samples, steps_between_samples, initial_n_chains;
	char sample_mode;
	config >> buff >> filename;
	config >> buff >> zmin;
	config >> buff >> zmax;
	config >> buff >> q;
	config >> buff >> target_swap_ratio;
	config >> buff >> swap_interval;
	config >> buff >> n_samples;
	config >> buff >> steps_between_samples;
	config >> buff >> initial_n_chains;
	config >> buff >> warmup;
	config >> buff >> sample_mode;
	config >> buff >> seed; 
	
	config.close();

	//set seed directly, then cal set_mmc for other parameters
	if (seed)
		srand(seed);
	else{
		seed = time(NULL);
	}
	set_mmc(zmin, zmax, q, target_swap_ratio, swap_interval, sample_mode, n_samples, steps_between_samples, initial_n_chains, warmup);
	add_initial_conformation_from_file(filename);
	return;
}

void mmchain::initialize(char* in, double zmin, double zmax, int q, double sr, int s, int n, int c, int nchains, int w, int m, int seed){
	string infile;
	infile.append(infile);
	set_mmc(zmin, zmax, q, sr, s, m, n, c, nchains, w);
	srand(seed);
	add_initial_conformation_from_file(infile);
}

void mmchain::create_config_file(){
	ofstream config("config.txt", ios::out);
	config << "initial_file= " << endl << "zmin= " << endl << "zmax= " << endl << "q= " << endl
		<< "target_swap_ratio= " << endl << "swap_interval= " << endl << "n_samples= " << endl 
		<< "steps_between_samples= " << endl << "initial_n_chains= 10" << endl
		<<"warmup= " << endl << "sample_mode= b"; 
	config.close();
}

bool mmchain::add_initial_conformation(istream& is){
	chain tempChain;
	int i=0;
	if (!initialComp0.readFromCoords(is))
      return false;
	 n_components = 1;
   if (initialComp1.readFromCoords(is))
      n_components = 2;

   if (n_components == 2){
	   for(i=0; i<m; i++){
		   tempChain.member_knot = new clkConformationBfacf3(initialComp0, initialComp1);
		   tempChain.member_knot->setSeed(time(NULL)); //what is this doing? PRESERVED(->setSeed(m+rand());)
		   chains.push_back(tempChain);
	   }}
   else if (n_components == 1){
		for (i=0;i<m;i++){
			tempChain.member_knot = new clkConformationBfacf3(initialComp0);
			tempChain.member_knot->setSeed(m+rand()); //what is this actually doing though?
			chains.push_back(tempChain);
		}}

   return true;
}

bool mmchain::add_initial_conformation_from_file(string& filename){ 
	ifstream in;
	in.open(filename.c_str(), ios::in);
	if (!in){
		cout << "ERROR: UNABLE TO OPEN FILE";
		return false;
	}
	if (add_initial_conformation(in))
		return true;
	else
		return false;
}

void mmchain::set_mmc(double Z_1, double Z_m, int Q, double Target_swap_ratio, int Swap_interval, char Sample_mode, int num_samples, int iterations, int M, double W){
	z_m = Z_m;
	z_1 = Z_1;
	q = Q;
	target_swap_ratio = Target_swap_ratio;
	swap_interval = Swap_interval;
	sample_mode = Sample_mode;
	n = num_samples;
	c = iterations;
	m = M;  
	w = W;
}

bool mmchain::stamp_log(string logname, string buff){
	logfile.open(logname.c_str(), std::ofstream::out | std::ofstream::app); //open file for appending
	if (logfile.fail()){
		cout << "Unable to create or open log.txt";
		return false;
	}
	logfile << buff << endl;
	//logfile << filename;
	logfile <<"zmin = "<< z_1 << endl;
	logfile <<"zmax= "<< z_m << endl;
	logfile <<"q= "<< q << endl;
	logfile <<"target_swap_ratio= " << target_swap_ratio << endl;
	logfile <<"swap_interval :"<<  swap_interval << endl;
	logfile <<"n_samples= " << n << endl;
	logfile <<"steps_between_samples :" << c << endl;
	logfile << "initial_n_chains: " << m << endl;
	logfile << "warmup= " << w << endl;
	logfile << "sample_mode= " << sample_mode << endl;
	logfile << "seed= " << seed << endl;
	cout << endl;
	logfile.close();

	return true;
}

void mmchain::run_mmc(){
	interval_data temp_interval;
	int i, j, k;
	double init_delta_b, b_1, b_m;
	//objects used only for testing lengths
	autocorr ac;
	autocorrInfo info;

	//add entry to log file
	std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);

    std::strftime(buffer,80,"%Y-%m-%d-%H-%M-%S",timeinfo);
    std::puts(buffer);
	stamp_log("log.txt", buffer);

	//open chronological block file
	char* buff = strcat(buffer, ".b");
	out.open(buff, std::ios::out);

	//calculate delta_b, produce initial spacing
	b_1 = log(z_1);
	b_m = log(z_m);
	init_delta_b = abs(b_1 - b_m)/m;

	//initialize chains and interval data
	chains[0].member_knot->setZ(z_m);
	for (i = 0; i < m; i++){
		chains[i].member_knot->setZ(exp((log(chains[0].member_knot->getZ()) - i*init_delta_b))); 
		chains[i].z = chains[i].member_knot->getZ();
		chains[i].data.resize(n);
	}
	for (i = 0; i < m-1; i++){
		temp_interval.delta_b = abs(log(chains[i].member_knot->getZ()) - log(chains[i+1].member_knot->getZ()));
		intervals.push_back(temp_interval); //consider static casting to save memory here
	}

	//warmup with swapping
	cout << "Warming up " << m << " chains: (" << w << " steps)" << endl;
	//update_size(); 
	for(int i = 0; i < w; i++){ 
		cout <<"\rCurrent Progress: " << i+1 << '/' << w;
		for(int j = 0; j < m; j++){
			chains[j].member_knot->stepQ(swap_interval, q, chains[j].z);
		}
		swap();
	}

	//calibrate chains
	calibrate_chains();
	cout << endl << "Starting sampling with " << m << " chains..." << endl;
	
	//main loop
	for(i=0; i<n; i++){
		for(j = 0; j < c/swap_interval; j++){
			for(k = 0; k < m; k++){
				chains[k].member_knot->stepQ(swap_interval, q, chains[k].z);
			}
			swap();
		}
		sample();
		cout << "\rSampling Progress: " << i+1 << '/' << n;
	}
	cout << endl;
	out.close();
	if ((sample_mode == 'b') || (sample_mode == 'a')){
		display_results();
	}
}

void mmchain::swap(){ 
	int i;
	clkConformationBfacf3* temp;
	probs* p_temp;
	if (n_components == 1){ //tested
		i = rand() % (m-1); //generate random integer from 0 to m-1
			if((chains[i].member_knot->getComponent(0).size() < chains[i+1].member_knot->getComponent(0).size()) 
				|| (((double) rand() / (RAND_MAX)) < pow((double)(chains[i+1].z)/(chains[i].z), (chains[i].member_knot->getComponent(0).size() - chains[i+1].member_knot->getComponent(0).size())))){
					//perform swap
				p_temp = chains[i].member_knot->probMap;
				chains[i].member_knot->probMap = chains[i+1].member_knot->probMap;
				chains[i+1].member_knot->probMap = p_temp;
				temp = chains[i].member_knot;
				chains[i].member_knot = chains[i+1].member_knot;
				chains[i+1].member_knot = temp;
				intervals[i].attempted_swaps++;
				intervals[i].sucessful_swaps++;
			}
			else{ 
				intervals[i].attempted_swaps++;
			}
		}
	if (n_components == 2){ //rewritten for new swap condition, untested
		for (i = 0; i < m-1;i++){
			i = rand() % (m-1); //generate random integer from 0 to m-1
			if((chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size()) 
				< (chains[i+1].member_knot->getComponent(0).size() + chains[i+1].member_knot->getComponent(1).size())
				|| (rand() < pow((double)(chains[i+1].z)/(chains[i].z), 
				((chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size()) 
				- (chains[i+1].member_knot->getComponent(0).size() + chains[i+1].member_knot->getComponent(1).size()))))){
				p_temp = chains[i].member_knot->probMap;
				chains[i].member_knot->probMap = chains[i+1].member_knot->probMap;
				chains[i+1].member_knot->probMap = p_temp;
					//perform swap
				temp = chains[i].member_knot;
				chains[i].member_knot = chains[i+1].member_knot;
				chains[i+1].member_knot = temp;
				intervals[i].attempted_swaps++;
				intervals[i].sucessful_swaps++;
			}
			else 
				intervals[i].attempted_swaps++;
		}
	}
}

chain* mmchain::alloc_new_chain(double z){
	// allocate and configure a single new chain from the file
	chain* tempChain = new chain();

	if (n_components == 2){
		tempChain->member_knot = new clkConformationBfacf3(initialComp0,  initialComp1);
	}
	else{
		tempChain->member_knot = new clkConformationBfacf3(initialComp0);
	}
	tempChain->member_knot->setZ(z);
	tempChain->z = tempChain->member_knot->getZ();
	//warmup new chain
	tempChain->member_knot->stepQ(w, q, tempChain->z); 
	return tempChain; 
}

void mmchain::insert_new_chain(int i){
	cout << "	Adding new chain at location: " << i+1 << endl;
	//allocate new chain/interval and adjust all other chain data according
	chain* tempChain = alloc_new_chain(exp((log(chains[i].z) + log(chains[i+1].z))/2));

	//insert new chain at position i+1, pushing all other chains back by one
	m += 1;
		chains.push_back(*tempChain); 
	for (int j = m-2; j > i; j--){
		chains[j+1] = chains[j];
	}
	chains[i+1] = *tempChain;
	//create new interval, push into invervals, update delta_bs using update_intervals()
	interval_data tempInterval;
	intervals.push_back(tempInterval); 
	for (int j = m-3; j > i; j--){
		intervals[j+1] = intervals[j];
	}
	intervals[i+1] = tempInterval;
	update_intervals(); 
}

void mmchain::update_intervals(){ //possibly obsolete
	for (int i = 0; i < m-1; i++){
		intervals[i].delta_b = (log(chains[i].z) - log(chains[i+1].z));
	}
}

bool mmchain::test_intervals(){
	int test = 1; //true by default
	for (int i = 0; i < m-1 ; i++){
		//compute confidence intervals
		double p = intervals[i].sucessful_swaps/intervals[i].attempted_swaps;
		double z = 1.96; //95% confidence
		double cw = z*sqrt((1/intervals[i].attempted_swaps)*p*(1-p));
		if ((p - cw) < target_swap_ratio){
			insert_new_chain(i);
			i++;
			test = 0; //false when new chain needs to be allocated
		}
		/*
		if (target_swap_ratio > p + cw){ //test confidence intervals
			//add new chain
			insert_new_chain(i);
			i++;
			test = 0; //false when new chain needs to be allocated
		}
		else if (target_swap_ratio < p - cw){
			system("PAUSE");
			//unhandled exception
		}*/
	}
	return test;
}

void mmchain::calibrate_chains(){
	cout << endl << "Starting calibration..." << endl;
	//Phase 1: Increase M as needed
	do{ 
		//reset interval data
		for (int i = 0; i < m-1; i++){
			intervals[i].attempted_swaps = 0;
			intervals[i].sucessful_swaps = 0;
		}
		//run chains for the same # of steps as a warmup
		cout << "\rTesting chains... " << endl;
		for(int i = 0; i < w; i++){ 
			for(int j = 0; j < m; j++){
				chains[j].member_knot->stepQ(swap_interval, q, chains[j].z);
			}
			swap();
		}
	}  while (!test_intervals()); //Compute and Check confidence intervals, add more chains as needed
}

void mmchain::sample(){
	if (n_components == 1){
		if(sample_mode == 's'){ //sample only
			for(int i = 0; i < m; i++){
				clkConformationAsList toPrint(chains[i].member_knot->getComponent(0));
				toPrint.writeAsCube(out);
			}
		}
		else if(sample_mode == 'a'){ //analyze length statistics
			for(int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size());
			}
		}
		else if(sample_mode == 'b'){ //both
			for(int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size());
				clkConformationAsList toPrint(chains[i].member_knot->getComponent(0));
				toPrint.writeAsCube(out);
			}
		}
	}
	else if (n_components == 2){
		if(sample_mode == 's'){ //sample only
			for(int i = 0; i < m; i++){
				clkConformationAsList toPrint(chains[i].member_knot->getComponent(0)), toPrint_2(chains[i].member_knot->getComponent(1));
				toPrint.writeAsCube(out);
				toPrint_2.writeAsCube(out);
			}
		}
		else if(sample_mode == 'a'){ //analyze length statistics
			for(int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size());
			}
		}
		else if(sample_mode == 'b'){ //both
			for(int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size());
				clkConformationAsList toPrint(chains[i].member_knot->getComponent(0)), toPrint_2(chains[i].member_knot->getComponent(1));
				toPrint.writeAsCube(out);
				toPrint_2.writeAsCube(out);
			}
		}
	}
}

void mmchain::display_results(){
	autocorr ac;
	autocorrInfo info;
	ofstream results;
	results.open("autocorr_data.txt", std::ofstream::out); //need to change the naming scheme to associate this log with the output file
	cout << endl;
	if ((sample_mode == 'b') || (sample_mode == 'a')){
	for (int i = 0; i < m; i++){
		info = ac.autocorrelation(chains[i].data,false);
		cout << chains[i].z << " " << info << endl;
		results << chains[i].z << " " << info << endl;	
	}}
	results.close();
}

void mmchain::update_size(){
		for (int i = 0; i < m; i++){
			chains[i].size = chains[i].member_knot->getComponent(0).size();
			cout << chains[i].size <<" ";
	}
		cout << endl;
}