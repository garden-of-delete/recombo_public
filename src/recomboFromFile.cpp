#include "recomboFromFile.h"

using namespace std;

recomboFromFile::recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int N_components, char Read_mode, int Sampling_mode, int Block_file_mode, bool Supress_output){
	//set operating variables
	supress_output = Supress_output;
	min_arc = Min_arc;
	max_arc = Max_arc;
	read_mode = Read_mode;
	sampling_mode = Sampling_mode;
	block_file_mode = Block_file_mode;
	//initialize infile_name and current_block_file_number, for use with block_file_mode == 1
	infile_name = new string(Infile);
	current_block_file_number = 0;
	//allocate uninitialized ifstream
	in = new ifstream;
	if (block_file_mode == 1){
		//incriment infile_name 
		inc_filename();
		//allocate and open output file
		out = new ofstream(Outfile, ios::out | ios::binary);
	}
	else{
		//open infile and allocate/open outfile
		in->open(Infile, ios::in | ios::binary);
		out = new ofstream(Outfile, ios::out | ios::binary);
	}
	n_components = N_components;
	if (!in->good())
		{cout << "Error reading input file" << endl; exit(3);}
	if (!out->good())
		{cout << "Error creating output file" << endl; exit(3);}
}

bool recomboFromFile::inc_filename(){
	//close current file if open
	in->close();
	//incriment block file counter
	current_block_file_number++;
	//put the new filename together with a string stream
	stringstream ss;
	ss << *infile_name << '%' << current_block_file_number << ".b";
	//try to open new filename
	in->open(ss.str().c_str(), ios::out | ios::binary);
	//check if succesful, return based on result
	if (in->good()){
		return true;
	}
	else{
		return false;
	}
}

string recomboFromFile::get_current_filename(){
	//put the new filename together with a string stream
	stringstream ss;
	ss << *infile_name << '%' << current_block_file_number << ".b";
	return ss.str();
}

bool recomboFromFile::read_comp_knots(ifstream* in){
	if (read_mode == 'b'){
		return initialComp0.readFromCube(*in);
	}
	else if (read_mode == 't'){
		return initialComp0.readFromText(*in);
	}
	else return false;
}

bool recomboFromFile::read_comp_links(ifstream* in){
	if (read_mode == 'b'){
		return (initialComp0.readFromCube(*in) && initialComp1.readFromCube(*in));
	}
	else if (read_mode == 't'){
		return (initialComp0.readFromText(*in) && initialComp1.readFromText(*in));
	}
	else return false;
}

void recomboFromFile::do_recombo(){
	if (sampling_mode == 0){
		sampling_mode = -1;
	}
	if (n_components == 1)
		do_recombo_knots();
	else if (n_components == 2)
		do_recombo_links();
	else{
	cout <<endl << "error: Links with > 2 components are not currently supported by this software" << endl;
	exit(2);
	}
}

void recomboFromFile::do_recombo_knots(){
	int attempts = 0, count = 0, sites = 0, choice = 0, length_counter = 0;
	long int total_attempts = 0, total_count = 0;
	vector<int> lengths;
	TOP:
	while (read_comp_knots(in) && sampling_mode != 0){
		sites = 0;
		if (!supress_output){
			cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
		}
		knot = new clkConformationBfacf3(initialComp0);
		//need to create new countRecomboSites function that is ideal for use with mmc
		int length = knot->getComponent(0).size();
		lengths.push_back(length);
		if (length == (min_arc + max_arc)){
		    length_counter ++;
            sites = knot->countRecomboSites(min_arc-1, max_arc+1);
        }
	//}
		//sites = knot->countRecomboSitesknot->getComponent(0).size()/2-4, knot->getComponent(0).size()/2+4);

		if(sites > 0){
			list<clkConformationAsList> components;
			choice = siteSelector.rand_integer(0, sites);
			knot->performRecombination();
			knot->getComponents(components);
			list<clkConformationAsList>::const_iterator i;
			for (i = components.begin(); i != components.end(); i++)
			{
				i->writeAsCube(*out);
			}
			count++;
			if (sampling_mode > 0){
				sampling_mode--;
			}
		}
		delete knot;
	}
	//report on current file
	cout << "\nPerformed "<<count<<" recombinations on "<<attempts<< " conformations in file "<< get_current_filename() << endl;
	double sum = 0;
	for (int i=0; i < lengths.size(); i++){
        sum += lengths[i];
	}
	cout << "Length " << min_arc + max_arc << " conformations so far: " << length_counter << endl;
	cout << "Average length so far: " << sum / lengths.size() << endl;
	if (block_file_mode == 1){
		//copy current file statistics into totals
		total_attempts += attempts;
		total_count += count;
		attempts = 0;
		count = 0;
		//if bfm == 1, see if next block file exists
		if (inc_filename()){
			//next block file exists and is now open, return to TOP to process next file
			goto TOP;
		}
		else{
			//next block file does not exist, report statistics for all blocks
			cout << "\nPerformed " << total_count << " recombinations on " << total_attempts << " conformations." << endl;
		}
	}
}


void recomboFromFile::do_recombo_links(){
	int attempts = 0, count = 0, sites = 0, choice = 0, length_counter = 0;
	long int total_attempts = 0, total_count = 0;
	vector<int> lengths;
	TOP:
	while (read_comp_links(in) && sampling_mode != 0){
	    sites = 0;
		if (!supress_output){
			cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
		}
		knot = new clkConformationBfacf3(initialComp0, initialComp1);
		//need to create new countRecomboSites function that is ideal for use with mmc
		int short_length = min(knot->getComponent(0).size(),knot->getComponent(1).size());
		int long_length = max(knot->getComponent(0).size(), knot->getComponent(1).size());
		lengths.push_back(short_length + long_length);
		if ((short_length >= min_arc) && (long_length <= max_arc)){
		    length_counter ++;
            sites = knot->countRecomboSites(min_arc-1, max_arc+1);
            }
	//}
		//sites = knot->countRecomboSites(knot->getComponent(0).size()/2-4, knot->getComponent(0).size()/2+4);

		if(sites > 0){
			list<clkConformationAsList> components;
			choice = siteSelector.rand_integer(0, sites);
			knot->performRecombination();
			knot->getComponents(components);
			list<clkConformationAsList>::const_iterator i;
			for (i = components.begin(); i != components.end(); i++)
			{
				i->writeAsCube(*out);
			}
			count++;
			if (sampling_mode > 0){
				sampling_mode--;
			}
		}
		delete knot;
	}
	//report on current file
	cout << "\nPerformed " << count << " recombinations on " << attempts << " conformations in file " << get_current_filename() << endl;
	double sum = 0;
	for (int i=0; i < lengths.size(); i++){
        sum += lengths[i];
	}
	cout << "Length " << min_arc + max_arc << " conformations so far: " << length_counter << endl;
	cout << "Average length so far: " << sum / lengths.size() << endl;
	if (block_file_mode == 1){
		//copy current file statistics into totals
		total_attempts += attempts;
		total_count += count;
		attempts = 0;
		count = 0;
		//if bfm == 1, see if next block file exists
		if (inc_filename()){
			//next block file exists and is now open, return to TOP to process next file
			goto TOP;
		}
		else{
			//next block file does not exist, report statistics for all blocks
			cout << "\nPerformed " << total_count << " recombinations on " << total_attempts << " conformations." << endl;
		}
	}
}

recomboFromFile::~recomboFromFile(){
	delete in;
	delete out;
	delete infile_name;
}
