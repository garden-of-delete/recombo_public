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

void mmchain::initialize(char* in, char* Outfile_name, double zmin, double zmax, int q, double sr, int s, int n, int c, long int w, int m, char mode, int seed, int bfs, bool Supress_output){
	string infile;
	min_arc = 0;
	max_arc = 0;
	target_recombo_length = 0;
	infile.append(in);
	set_mmc(zmin, zmax, q, sr, s, mode, n, c, m, w);
	srand(seed);
	add_initial_conformation_from_file(infile);
	block_file_size = bfs;
	block_file_index = 0;
	outfile_name = Outfile_name;
	current_block_file_number = 0;
	supress_output = Supress_output;
	sample_attempts = 0;
}

void mmchain::initialize(char* in, char* Outfile_name, double zmin, double zmax, int q, double sr, int s, int n, int c, long int w, int m, char mode, int Seed,
	int Min_arc, int Max_arc, int Target_recombo_length, int bfs, bool Supress_output){
	if (Min_arc == 0 && Max_arc == 0 && Target_recombo_length == 0){ //if not in filtering mode, use simpler constructor
		initialize(in, Outfile_name, zmin, zmax, q, sr, s, n, c, w, m, mode, Seed, bfs, Supress_output);
		return;
	}
	string infile;
	min_arc = Min_arc;
	max_arc = Max_arc;
	target_recombo_length = Target_recombo_length;
	infile.append(in);
	set_mmc(zmin, zmax, q, sr, s, mode, n, c, m, w);
	seed = Seed;
	srand(seed);
	add_initial_conformation_from_file(infile);
	block_file_size = bfs;
	block_file_index = 0;
	outfile_name = Outfile_name;
	current_block_file_number = 0;
	supress_output = Supress_output;
	sample_attempts = 0;
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

void mmchain::set_mmc(double Z_1, double Z_m, int Q, double Target_swap_ratio, int Swap_interval, char Sample_mode, int num_samples, int iterations, int M, long int W){
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

void mmchain::stamp_log(string buff){
	
	cout << "start_time= " << buff << endl;
	cout <<"zmin= "<< z_1 << endl;
	cout <<"zmax= "<< z_m << endl;
	cout <<"q= "<< q << endl;
	cout <<"target_swap_ratio= " << target_swap_ratio << endl;
	cout <<"swap_interval :"<<  swap_interval << endl;
	cout <<"n_samples= " << n << endl;
	cout <<"steps_between_samples :" << c << endl;
	cout << "initial_n_chains: " << m << endl;
	cout << "warmup= " << w << endl;
	cout << "sample_mode= " << sample_mode << endl;
	cout << "seed= " << seed << endl;
	cout << "block_file_size= " << block_file_size << endl;
	cout << endl;
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
    //std::puts(buffer); moved inside stamp_log(...)
	stamp_log(buffer);

	/*
	//open chronological block file (legacy code, here for reference)
	char* buff = strcat(buffer, ".b");
	out.open(buff, std::ios::out | std::iostream::binary);
	*/
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
		intervals.push_back(temp_interval);
	}

	//warmup with swapping
	cout << "Warming up " << m << " chains: (" << w << " steps)" << endl;
	//update_size(); 
	for(int i = 0; i < w; i++){ 
		if (!supress_output){
			cout << "\rCurrent Progress: " << i + 1 << '/' << w;
		}
		for(int j = 0; j < m; j++){
			chains[j].member_knot->stepQ(swap_interval, q, chains[j].z);
		}
		swap();
	}

	//calibrate chains
	calibrate_chains();
	cout << endl << "Starting sampling with " << m << " chains..." << endl;
	
	//main loop
	i = 0;
	while(i < n || sample_mode == 'e'){
		cout << i << endl;
		for(j = 0; j < c/swap_interval; j++){
			for(k = 0; k < m; k++){
				chains[k].member_knot->stepQ(swap_interval, q, chains[k].z);
			}
			swap();
		}
		i += sample();
		sample_attempts++;
		if (sample_mode != 'f'){
			i++;
		}
		if (!supress_output){
			cout << "\rSampling Progress: " << i << '/' << n;
		}
	}
	cout << endl;
	if (sample_mode == 'f'){
		cout << "Sampled " << i << " times from " << sample_attempts << " attempts." << endl;
		cout << "Ratio: " << float(i) / float(sample_attempts) << endl;
	}
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
			//test code
			//int shouldBeLarger = chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size();
			//int shouldBeSmaller = chains[i + 1].member_knot->getComponent(0).size() + chains[i + 1].member_knot->getComponent(1).size();
			if((chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size()) 
				< (chains[i+1].member_knot->getComponent(0).size() + chains[i+1].member_knot->getComponent(1).size())
				|| (((double)rand() / (RAND_MAX)) < pow((double)(chains[i + 1].z) / (chains[i].z),
				((chains[i].member_knot->getComponent(0).size() + chains[i].member_knot->getComponent(1).size()) 
				- (chains[i+1].member_knot->getComponent(0).size() + chains[i+1].member_knot->getComponent(1).size()))))){
				//swap pre-computed bfacf probabilities
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
		cout << "Interval: " << i << ": " << p - cw << ' ';
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

int mmchain::sample(){
	if (n_components == 1){
		if (sample_mode == 's'){ //sample only
			for (int i = 0; i < m; i++){
				write_to_block_file(chains[i].member_knot);
			}
		}
		else if (sample_mode == 'a'){ //analyze length statistics
			for (int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size());
			}
		}
		else if (sample_mode == 'b'){ //both
			for (int i = 0; i < m; i++){
				chains[i].data.push_back(chains[i].member_knot->getComponent(0).size());
				write_to_block_file(chains[i].member_knot);
			}
		}
		else if (sample_mode == 'f'){ //filter sampling
			int samples = 0;
			for (int i = 0; i < m; i++){
				if (do_recombo_knots(i))
					samples++;
				}
			return samples;
			}
		else if (sample_mode == 'e'){ //eternal sampling mode
			for (int i = 0; i < m; i++){
				write_to_block_file(chains[i].member_knot);
			}
		}
	}
	else if (n_components == 2){
		if(sample_mode == 's'){ //sample only
			for(int i = 0; i < m; i++){
				clkConformationAsList toPrint(chains[i].member_knot->getComponent(0)), toPrint_2(chains[i].member_knot->getComponent(1));
				write_to_block_file(chains[i].member_knot);
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
				write_to_block_file(chains[i].member_knot);
			}
		}
		else if (sample_mode == 'f'){ //filter sampling
			int samples = 0;
			for (int i = 0; i < m; i++){
				if (do_recombo_links(i))
					samples++;
			}
			return samples;
		}
		else if (sample_mode == 'e'){ //eternal sampling mode
			for (int i = 0; i < m; i++){
				write_to_block_file(chains[i].member_knot);
			}
		}
	}
	return 0;
}

bool mmchain::do_recombo_knots(int current_chain){
	//setup recombo environment
	int sites = 0, choice = 0;
	pseudorandom siteSelector;

	//perform recombination
		//check length and count recombo sites
	int length = chains[current_chain].member_knot->getComponent(0).size();
	if (length == (min_arc + max_arc)){
		sites = chains[current_chain].member_knot->countRecomboSites(min_arc - 1, max_arc + 1);
	}
	if (sites > 0){
		choice = siteSelector.rand_integer(0, sites-1);
		//chains[current_chain].member_knot->performRecombination(choice);
			//perform recombination and write both conformations to respective block files
		write_recombination_to_block_file(chains[current_chain].member_knot, choice);
		//chains[current_chain].member_knot->undoRecombination();
		return true;
	}
	return false;
}

bool mmchain::do_recombo_links(int current_chain){ 
	//setup recombo environment
	int sites = 0, choice = 0;
	pseudorandom siteSelector;

	//perform recombination
	//check length and count recombo sites
	int length = chains[current_chain].member_knot->getComponent(0).size();
	if (length == (min_arc + max_arc)){
		sites = chains[current_chain].member_knot->countRecomboSites(min_arc - 1, max_arc + 1);
	}
	if (sites > 0){
		choice = siteSelector.rand_integer(0, sites - 1);
		write_recombination_to_block_file(chains[current_chain].member_knot, choice);
		return true;
	}
	return false;
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

void mmchain::write_to_block_file(clkConformationBfacf3* clk){
	if ((block_file_index < block_file_size) && (block_file_index != 0)){
		//write conformation to file
		if (n_components == 1){
			clkConformationAsList toPrint(clk->getComponent(0));
			toPrint.writeAsCube(out);
		}
		else if (n_components == 2){
			clkConformationAsList toPrint(clk->getComponent(0));
			clkConformationAsList toPrint2(clk->getComponent(1));
			toPrint.writeAsCube(out);
			toPrint2.writeAsCube(out);
		}
		else{
			cout << endl << "ERROR: write_to_block_file(): n_components=" << n_components << ", only 1 or 2 component links supported" << endl;
			exit(1);
		}
		block_file_index++;
	}
	else 
	{
		//close current file if open
		out.close();
		//incriment file name
		current_block_file_number++;
		//reset block file index
		block_file_index = 0;
		stringstream ss;
		ss << outfile_name << '%' << current_block_file_number << ".b";
		//open new file
		out.open(ss.str().c_str(), ios::out | ios::binary);
		
		//write conformation to file
		if (n_components == 1){
			clkConformationAsList toPrint(clk->getComponent(0));
			toPrint.writeAsCube(out);
		}
		else if (n_components == 2){
			clkConformationAsList toPrint(clk->getComponent(0));
			clkConformationAsList toPrint2(clk->getComponent(1));
			toPrint.writeAsCube(out);
			toPrint2.writeAsCube(out);
		}
		else{
			cout << endl << "ERROR: write_to_block_file(): n_components=" << n_components << ", only 1 or 2 component links supported" << endl;
			exit(1);
		}
		block_file_index++;
	}
}

void mmchain::write_recombination_to_block_file(clkConformationBfacf3* clk, int site_choice){
	if ((block_file_index < block_file_size) && (block_file_index != 0)){
		if (n_components == 1){
			//write before recombination conformation to file
			clkConformationAsList toPrint(clk->getComponent(0));
			toPrint.writeAsCube(out);
			//perform recombination
			clk->performRecombination(site_choice);
			//write clk_after to second output file
			clkConformationAsList toPrint2(clk->getComponent(0));
			toPrint.writeAsCube(out2);
			//undo recombination
			clk->undoRecombination();
		}
		else if (n_components == 2){
			//write before recombination conformation to file
			clkConformationAsList toPrint(clk->getComponent(0));
			clkConformationAsList toPrint2(clk->getComponent(1));
			toPrint.writeAsCube(out);
			toPrint2.writeAsCube(out);
			//perform recombination
			clk->performRecombination(site_choice);
			//write clk_after to second output file
			clkConformationAsList toPrint3(clk->getComponent(0));
			clkConformationAsList toPrint4(clk->getComponent(1));
			toPrint3.writeAsCube(out2);
			toPrint4.writeAsCube(out2);
			//undo recombination
			clk->undoRecombination();
		}
		else{
			cout << endl << "ERROR: write_to_block_file(): n_components=" << n_components << ", only 1 or 2 component links supported" << endl;
			exit(1);
		}
		block_file_index++;
	}
	else
	{
		//close current file if open
		out.close();
		//incriment file name
		current_block_file_number++;
		//reset block file index
		block_file_index = 0;
		stringstream ss;
		ss << outfile_name << '%' << current_block_file_number << ".b";
		//open new file
		out.open(ss.str().c_str(), ios::out | ios::binary);

		//if in filter mode, open second output file
		if (sample_mode == 'f'){ //this conditional block is probably not neccisary 
			stringstream tt;
			tt << outfile_name << "_after%" << current_block_file_number << ".b";
			out2.open(tt.str().c_str(), ios::out | ios::binary);
		}

		if (n_components == 1){
			//write before recombination conformation to file
			clkConformationAsList toPrint(clk->getComponent(0));
			toPrint.writeAsCube(out);
			//perform recombination
			clk->performRecombination(site_choice);
			//write clk_after to second output file
			clkConformationAsList toPrint2(clk->getComponent(0));
			toPrint.writeAsCube(out2);
			//undo recombination
			clk->undoRecombination();
		}
		else if (n_components == 2){
			//write before recombination conformation to file
			clkConformationAsList toPrint(clk->getComponent(0));
			clkConformationAsList toPrint2(clk->getComponent(1));
			toPrint.writeAsCube(out);
			toPrint2.writeAsCube(out);
			//perform recombination
			clk->performRecombination(site_choice);
			//write clk_after to second output file
			clkConformationAsList toPrint3(clk->getComponent(0));
			clkConformationAsList toPrint4(clk->getComponent(1));
			toPrint3.writeAsCube(out2);
			toPrint4.writeAsCube(out2);
			//undo recombination
			clk->undoRecombination();
		}
		else{
			cout << endl << "ERROR: write_to_block_file(): n_components=" << n_components << ", only 1 or 2 component links supported" << endl;
			exit(1);
		}
		block_file_index++;
	}
}
