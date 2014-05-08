#include "recomboFromFile.h"


recomboFromFile::recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int N_components/*, char Mode*/){
	min_arc = Min_arc;
	max_arc = Max_arc;
	in = new ifstream(Infile, ios::in | ios::binary);
	out = new ofstream(Outfile, ios::out | ios::binary);
	//mode = Mode;
	n_components = N_components;
	if (!in->good())
		{cout << "Error reading input file" << endl; exit(3);}
	if (!out->good())
		{cout << "Error reading output file" << endl; exit(3);}

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
	while (initialComp0.readFromCube(*in)){
	    sites = 0;
		cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
		knot = new clkConformationBfacf3(initialComp0);
		//need to create new countRecomboSites function that is ideal for use with mmc
		int length = knot->getComponent(0).size();
		lengths.push_back(length);
		if (length == 100){
		    length_counter ++;
            sites = knot->countRecomboSites(48, 52);
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
	cout << "Length 100 conformations: " << length_counter << endl;
	cout << "Average Length: " << sum / lengths.size() << endl;
}


void recomboFromFile::do_recombo_links(){
	int attempts = 0, count = 0, sites = 0, choice = 0, length_counter = 0;
	vector<int> lengths;
	while (initialComp0.readFromCube(*in) && initialComp1.readFromCube(*in)){
	    sites = 0;
		cout << '\r' << "Attempting Recombination: " << ++attempts << " Performed: " << count;
		knot = new clkConformationBfacf3(initialComp0, initialComp1);
		//need to create new countRecomboSites function that is ideal for use with mmc
		int length = knot->getComponent(0).size() + knot->getComponent(1).size();
		lengths.push_back(length);
		if (length == 100){
		    length_counter ++;
			knot->
            sites = knot->countRecomboSites(48, 52);
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
	cout << "Length 100 conformations: " << length_counter << endl;
	cout << "Average Length: " << sum / lengths.size() << endl;
}

recomboFromFile::~recomboFromFile(){
	delete in;
	delete out;
}