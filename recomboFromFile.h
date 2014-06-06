#include "mmchain.h"

class recomboFromFile{
private:
	int min_arc, max_arc, n_components;
	clkConformationAsList initialComp0, initialComp1;
	clkConformationBfacf3* knot;
	pseudorandom siteSelector;
	ifstream* in;
	ofstream* out;
	void do_recombo_knots();
	void do_recombo_links();
public:
	recomboFromFile(int Min_arc, int Max_arc, char* Infile, char* Outfile, int n_components);
	void do_recombo();
	~recomboFromFile();
};

