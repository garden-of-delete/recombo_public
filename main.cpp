#include <iostream>
#include <list>
#include <fstream>

//#include "mmchain.h"
#include "random.h"
#include "autocorr.h"
#include "clk.h"
#include "clkCigar.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "conformationAsList.h"
#include "genericConformation.h"
#include "Analyzer.h"
//#include "mmchain.h" included in recomboFromFile.h
#include <string.h>
#include "runRecombo.h"
#include "recomboArgprocessor.h"
#include "recomboFromFile.h"


#include "argprocessor.h"
//#include "asdf.h"

using namespace std;
/*
int main(void){
    mmchain mmc;
    mmc.initialize();
    mmc.run_mmc();
}*/
/*
int main(void){
	
    recomboFromFile recombo(50, 50, "4_2_1a_before.b", "4_2_1a_after_strict.b", 2);
    recombo.do_recombo();
    system("PAUSE");
}*/

int main(void){
	clkConformationBfacf3* knot;
	vector<double> lengths;
	ifstream is;
	is.open("initial/unknot");
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	knot = new clkConformationBfacf3(initialComp0);
	double z = .2;
	int q = 3;
	knot->init_Q(.2000,3);
	int p = 1000000000;
	for (int i=0; i < p; i++){
		knot->stepQ(q,z);
		cout << '\r' << i+1 << '/' << p;
	}

}
/*
//Q Test code
int main(void){
	clkConformationBfacf3* knot;
	vector<double> data;
	ifstream is;
	is.open("initial/unknot");
	//clkConformationAsList convert;
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	double means[20], autocorrs[20], mean_vars[20], autocorr_vars[20];
    	double z = .2175;
	srand(478);
	for (int i=0; i < 1; i++){
        knot = new clkConformationBfacf3(initialComp0);
        int q = 3;
        //knot->init_Q(knot->getComponent(0).size(), knot->getZ(),3,false);
        knot->stepQ(5000000000,q,z);
	cout << "Warmup Complete!" << endl;
        for (int i = 0; i < 10000; i++){
            knot->stepQ(1000000,q,z);
            data.push_back(knot->getComponent(0).size());
            //cout << data.size() <<' ';
        }
            autocorr ac;
            autocorrInfo info;
        info = ac.autocorrelation(data, false);
        cout << "With Q = " << q << ":" << endl << info << endl;
        means[i] = info.mean;
        mean_vars[i] = info.error_mean;
        autocorrs[i] = info.autocorr;
        autocorr_vars[i] = info.error_autocorr;
        data.clear();
        delete knot;
	}
	double sum0=0, sum1=0, sum2=0, sum3=0;
	for(int i=0; i < 20; i++){
        sum2 += means[i];
        sum3 += mean_vars[i];
        sum0 += autocorrs[i];
        sum1 += autocorr_vars[i];
	}
	sum0 /= 20;
	sum1 /= 20;
	sum2 /= 20;
	sum3 /= 20;

	cout <<"z:" << z << " Autocorr Average: " << sum0 << " Autocorr stdev: " << sqrt(sum1) << " Mean Length: " << sum2 << "Length stdev: " << sqrt(sum3) << endl;
}*/
