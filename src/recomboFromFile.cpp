#include "recomboFromFile.h"

using namespace std;

recomboFromFile::recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int N_components, char Read_mode){
	min_arc = Min_arc;
	max_arc = Max_arc;
	read_mode = Read_mode;
	in = new ifstream(Infile, ios::in | ios::binary);
	out = new ofstream(Outfile, ios::out | ios::binary);
	//mode = Mode;
	n_components = N_components;
	if (!in->good())
		{cout << "Error reading input file" << endl; exit(3);}
	if (!out->good())
		{cout << "Error creating output file" << endl; exit(3);}
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
	vector<int> lengths;
	while (read_comp_knots(in)){
	    sites = 0;
		cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
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
		}
		delete knot;
	}
	cout << "\nPerformed "<<count<<" recombinations on "<<attempts << " conformations."<< endl;
	double sum = 0;
	for (int i=0; i < lengths.size(); i++){
        sum += lengths[i];
	}
	cout << "Length 100 conformations: " << length_counter << endl;
	cout << "Average Length: " << sum / lengths.size() << endl;
}


void recomboFromFile::do_recombo_links(){
	int attempts = 0, count = 0, sites = 0, choice = 0, length_counter = 0;
	vector<int> lengths;
	while (read_comp_links(in)){
	    sites = 0;
		cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
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
		}
		delete knot;
	}
	cout << "\nPerformed "<<count<<" recombinations on "<<attempts << " conformations."<< endl;
	double sum = 0;
	for (int i=0; i < lengths.size(); i++){
        sum += lengths[i];
	}
	cout << "Length 50/50 conformations: " << length_counter << endl;
	cout << "Average Length: " << sum / lengths.size() << endl;
}

recomboFromFile::~recomboFromFile(){
	delete in;
	delete out;
}
