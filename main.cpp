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
#include "Analyzer.h"
//#include "mmchain.h" included in recomboFromFile.h
#include <string.h>
#include "runRecombo.h"
#include "recomboArgprocessor.h"
#include "recomboFromFile.h"


#include "argprocessor.h"

//#include "asdf.h"

using namespace std;

bool autocorrLength(const clk& initial, double z, int seed, int& c, double& meanlength)
{
   vector<double> data;
   autocorr ac;
   autocorrInfo info;
   info.autocorr = info.error_autocorr = 10.0;
   clkConformationBfacf3 knot(initial);
   knot.setZ(z);
   knot.setSeed(seed);
   int timeout = 20; // try loop this number of times at most
   // make initial guess for c and samplesize
   c = 5000;
   int samplesize = 2000;

   while (timeout && (info.error_autocorr > 0.01 || info.autocorr > 0.55 || info.error_mean / abs(info.mean) > 0.1))
   {
      timeout--;
      data.clear();
      cout << "Taking " << samplesize << " samples every " << c << " iterations." << endl;
      for (int i = 0; i < samplesize; i++)
      {
         knot.step(c);
         data.push_back(knot.getComponent(0).size());
      }
      info = ac.autocorrelation(data, false);
      //      cout << info << endl;
      if (info.autocorr > 0.6)
      {
         int d = 2 * c * info.autocorr;
         if (d > 2 * c) c *= 2;
         else c = d;
      }
      if (info.error_autocorr > 0.01 || info.error_mean / abs(info.mean) > 0.05)
      {
         samplesize *= 2;
         cout << "Increasing samplesize to " << samplesize << "." << endl;
      }
   }

   return timeout != 0;
}

void searchForLongConformations(const clk& initial, double z, int targetLength, int seed, int c, int desired)
{
   //   vector<double> greaterLength, exactLength;
   //   autocorr ac;
   //   autocorrInfo info;
   ofstream out("l100.b", ios::binary | ios::out);
   clkConformationBfacf3 knot(initial);
   knot.setZ(z);
   knot.setSeed(seed);

   pseudorandom s(seed * 923589);

   // do a warmup
   knot.step(20 * c);

   int i = 0;
   while (desired)
   {
      // take no samples from c iterations
      knot.step(2 * c);

      // Now save the state of random number generator
      pseudorandom p(saveRandomState());
      // make two copies of conformation so we can duplicate a run
      clkConformationAsList savedConformation(knot.getComponent(0));
      threevector<int> v;
      savedConformation.getVertex(0, v);
      v *= -1;
      savedConformation.translate(v);
      //      cout << savedConformation;
      //      clkConformationBfacf3 knot1(knot.getComponent(0)), knot2(knot.getComponent(0));
      clkConformationBfacf3 knot1(savedConformation);

      // now search for conformations of requested length or greater than requested
      int greater, exact;
      greater = exact = 0;
      for (int j = 0; j < c; j++)
      {
         knot1.step();
         if (knot1.getComponent(0).size() == targetLength)
            exact++;
         if (knot1.getComponent(0).size() >= targetLength)
            greater++;
      }
      if (exact == 0) cout << i << " Nope, nuttin'." << endl;
      else
      {
         // if we are here, then we have found at least one conformation matching
         // criterion.
         desired--;
         bool found = false; // this is for a paranoia check

         cout << i << " exact = " << exact << ", and greater = " << greater
                 << "." << endl;
         // Save the state of the random number generator so that we can restore later
         pseudorandom q(saveRandomState());

         // restore random number generator to redo run
         copyRandomState(p);

         // decide which of the conformations we will choose
         int choice = s.rand_integer(1, exact);

         clkConformationBfacf3 knot2(savedConformation);
         for (int j = 0; j < c; j++)
         {
            knot2.step();
            if (knot2.getComponent(0).size() == targetLength)
            {
               choice--;
               if (choice == 0)
               {
                  // we have the conformation we wanted
                  found = true;
                  clkConformationAsList chosenConformation(knot2.getComponent(0));
                  chosenConformation.getVertex(0, v);
                  v *= -1;
                  chosenConformation.translate(v);
                  cout << chosenConformation.size() << endl;
                  cout << chosenConformation << endl;
                  chosenConformation.writeAsCube(out);
                  break;
               }
            }
         }

         if (!found)
         {
            cout << "Oh, no." << endl;
            return;
         }

         // restore random number generator again
         copyRandomState(q);
      }

      i++;
   }

   out.close();
}

// TODO: move this to unit tests

void hong()
{
   double w1 = 137.27;
   double w2 = 138.643;
   double w3 = 140.595;

   double ch1[3 * 32] = {-37.3989606529101, -3.10588828391642, 20.586042721413, -37.4717276074902, -3.6396136407412, 19.7435213781799, -38.0202408491441, -3.9271891370853, 20.5286542061987, -37.9580233910169, -3.79701952376605, 19.5391165062233, -37.6706787138507, -3.12933642492661, 18.8523615892199, -37.365006464846, -2.36307983983314, 18.2871931877772, -36.8013940399895, -1.57085051203344, 18.0532835327362, -36.6659351714447, -1.37629989812607, 19.0247777188025, -36.4299802291331, -1.08367118753633, 19.9514349925221, -36.4492736062506, -0.114048552870768, 19.7075913706171, -36.6930683869673, -0.616283716449046, 20.5372442771755, -36.7661557890101, -0.77550507777981, 21.5217780702767, -37.2328686756442, -1.65977185786103, 21.5376318740662, -36.4249647519895, -1.15451757770988, 21.2343014230192, -36.5822211747878, -1.33646193058231, 22.2049541260735, -36.9849670532546, -1.77111008881778, 21.3994252712705, -37.5807363960653, -2.39261264527729, 20.8907078030414, -38.2329243966089, -3.13206337610574, 21.0576312244787, -37.3784492682881, -2.96066859234085, 20.5672272245159, -37.7980133266929, -2.19787854151305, 20.0751730209436, -36.9154043517058, -1.91662115642744, 20.4518634516783, -37.6364997098321, -2.16590194165702, 19.8054266221241, -36.7903828211506, -2.69579579478075, 19.8628610639675, -37.6392091845022, -2.41690571801364, 19.4137351017249, -37.6221789713701, -2.64743361681578, 18.4408184250692, -38.050021186594, -2.28166459021156, 17.6142810344756, -37.1877871502802, -1.78331915559176, 17.7048585946981, -36.7590232994614, -1.59404093680091, 18.5882244963055, -37.6799876180177, -1.96381368690736, 18.4653715874068, -36.8999618864412, -1.74431522034427, 19.0513581753549, -37.7724924358349, -1.61529047974095, 19.5225726240796, -37.2805835136952, -2.41037568121079, 19.8773473902665};

   double ch2[3 * 24] = {619.397625517108, 28.1668454960761, 179.256207098995, 620.219596717076, 27.7041728551382, 179.588317583291, 621.028887401871, 27.3551451506494, 180.060787951679, 621.94999935394, 26.9744766773815, 179.979275466758, 622.912253832248, 26.7274537986462, 179.865056307416, 623.902094985486, 26.8519585928097, 179.796404667959, 624.682221962914, 27.4769221462741, 179.82508316003, 624.967279289255, 28.1039889252456, 179.100148973088, 625.476770214304, 28.9524009385221, 178.956566059655, 625.527626048141, 29.6786187599647, 178.270985021462, 624.823114997854, 30.1233873702694, 177.717952415426, 624.27621388262, 29.9545606208765, 176.897954428311, 623.648033696763, 29.6632348700038, 176.176484428237, 623.070205460253, 29.2812994389729, 175.455207882372, 622.387176600236, 29.3437100482762, 174.727487817897, 621.515961391931, 29.1353692318457, 174.282990163162, 620.615340693876, 29.1751975163691, 173.850213008002, 619.657365215002, 29.2963867333111, 174.110205614772, 618.825751052682, 29.2747019203403, 174.665135929273, 618.124178426531, 28.8635295325945, 175.247143762366, 617.633125779214, 28.4808347742654, 176.029711343819, 617.66861258145, 28.3226271969996, 177.016479325343, 618.36694149694, 28.1403865320039, 177.708667930996, 619.097768757941, 28.4128562662067, 178.334488798764};

   double ch3[3 * 100] = {20, 5, -14, 20, 5, -15, 21, 5, -15, 21, 5, -16, 21, 4, -16, 21, 3, -16, 22, 3, -16, 22, 4, -16, 22, 4, -15, 23, 4, -15, 24, 4, -15, 25, 4, -15, 25, 4, -16, 25, 5, -16, 25, 5, -15, 25, 6, -15, 24, 6, -15, 23, 6, -15, 23, 5, -15, 23, 5, -14, 22, 5, -14, 22, 6, -14, 22, 6, -13, 22, 5, -13, 21, 5, -13, 21, 5, -12, 20, 5, -12, 20, 6, -12, 20, 6, -13, 20, 5, -13, 20, 4, -13, 20, 4, -14, 20, 3, -14, 20, 2, -14, 19, 2, -14, 19, 3, -14, 19, 4, -14, 18, 4, -14, 18, 4, -15, 18, 3, -15, 19, 3, -15, 19, 4, -15, 19, 5, -15, 18, 5, -15, 18, 5, -16, 17, 5, -16, 17, 6, -16, 18, 6, -16, 19, 6, -16, 20, 6, -16, 20, 5, -16, 20, 5, -17, 20, 5, -18, 20, 5, -19, 20, 6, -19, 21, 6, -19, 22, 6, -19, 22, 5, -19, 23, 5, -19, 24, 5, -19, 24, 5, -18, 24, 4, -18, 24, 4, -17, 24, 4, -16, 24, 5, -16, 24, 5, -15, 24, 5, -14, 24, 4, -14, 24, 4, -13, 24, 5, -13, 24, 6, -13, 23, 6, -13, 23, 6, -12, 23, 7, -12, 23, 8, -12, 23, 9, -12, 23, 9, -11, 22, 9, -11, 21, 9, -11, 21, 8, -11, 21, 8, -10, 20, 8, -10, 20, 8, -11, 20, 8, -12, 20, 8, -13, 19, 8, -13, 19, 8, -12, 19, 7, -12, 19, 7, -13, 20, 7, -13, 21, 7, -13, 21, 7, -14, 21, 7, -15, 20, 7, -15, 20, 8, -15, 19, 8, -15, 19, 7, -15, 19, 6, -15, 19, 6, -14, 19, 5, -14};

   conformationAsList con1, con2, con3;



   for (int i = 0; i < 32; i++)
   {
      con1.addVertexBack(ch1[3 * i], ch1[3 * i + 1], ch1[3 * i + 2]);
   }

   for (int i = 0; i < 24; i++)
   {
      con2.addVertexBack(ch2[3 * i], ch2[3 * i + 1], ch2[3 * i + 2]);
   }

   for (int i = 0; i < 100; i++)
   {
      con3.addVertexBack(ch3[3 * i], ch3[3 * i + 1], ch3[3 * i + 2]);
   }

   //   cout << con1 << endl << endl << con2 << endl << endl << con3 << endl;

   double Wr, Acn;
   con1.writheACN(Wr, Acn);
   //   computeWrAcn(con1.data.begin(), con1.data.end(), Wr, Acn);
   cout << Wr << endl;
   //   cout << con1.rog() << endl;

   con2.writheACN(Wr, Acn);
   //   computeWrAcn(con2.data.begin(), con2.data.end(), Wr, Acn);
   cout << Wr << endl;

   con3.writheACN(Wr, Acn);
   //   computeWrAcn(con3.data.begin(), con3.data.end(), Wr, Acn);
   cout << Wr << endl;
}

//int mainHong();

int mainRecombo(int argc, char** argv);
int mainRerecombo(int argc, char** argv);
int mainTwoCompWrithe(int argc, char** argv);
int mainFixEgc(int argc, char** argv);

int main(int argc, char**argv)
{
	//srand(time(NULL));
	//Analyzer Test Code
	//Analyzer analyzer;
	//analyzer.add_initial_conformation_from_file("coords/3_1");
	//analyzer.z_from_length(100, .5);
	//analyzer.generate_table("initial/3_1", "results/3_1.csv", 50, 70, 1);


//mmc test code
	
  mmchain test_chain;
  test_chain.initialize();
  test_chain.run_mmc();
  //test_chain.display_results();*/

	
//recomboFromFile test code
	/*recomboFromFile recombo(50, 50, "3_1_before.b", "3_1_after.b", 1);
	recombo.do_recombo();*/
	
		/*
//clkConformationAsList.writeAsCube() test
	clkConformationAsList test;
	ifstream coords("3_1", ios::in);
	ofstream out("3_1_cube", ios::out);
	test.readFromCoords(coords);
	test.writeAsCube(out);
	out.close();
	test.clear();
	ifstream in("3_1_cube", ios::in | ios::binary);
	test.readFromCube(in);
	in.close();*/


  /*//comparison test code
	clkConformationBfacf3* knot;
	vector<double> data;
	ifstream is;
	is.open("5_1");
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	knot = new clkConformationBfacf3(initialComp0);
	knot->getComponent(0).setZ(.19273);
	knot->setSeed(42);
	clkConformationAsList one(knot->getComponent(0));
	for (int i = 0; i < 5000; i++){
		knot->step(10000);
		data.push_back(knot->getComponent(0).size());
	}
	autocorr ac;
	autocorrInfo info;
	info = ac.autocorrelation(data, false);
	cout << endl <<"Without MMC:"<<endl;
	cout << info;
	*/

	/*
	//Sampling code for reika
	int j=1;
    clkConformationBfacf3* knot;
	vector<double> data;
	ifstream is;
	is.open("initial/7_5");
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	knot = new clkConformationBfacf3(initialComp0);
	knot->getComponent(0).setZ(.18981);
	knot->setSeed(41);
	string filename;
	while (j <= 100){
		knot->step(250000); 
		if (knot->getComponent(0).size() == 100){
			stringstream convert;
			ofstream out;
			convert << j;
			filename = convert.str();
			out.open("results/samples/7_5/"+filename);
			clkConformationAsList toPrint(knot->getComponent(0));
			toPrint.writeAsText(out);
			out.close();
			cout << " " << j << endl;
			j++;
		}
	}
	return 0;*/

	/*
	//Q Test code
	clkConformationBfacf3* knot; 
	vector<double> data;
	ifstream is;
	is.open("initial/unknot");
	//clkConformationAsList convert;
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	knot = new clkConformationBfacf3(initialComp0);
	double z = .2100;
	int q = 3;
	knot->setSeed(54); //set seed for test purposes
	//knot->init_Q(knot->getComponent(0).size(), knot->getZ(),3,false);
	knot->stepQ(1000000000,q,z);
	for (int i = 0; i < 5000; i++){
		knot->stepQ(900000,q,z);
		data.push_back(knot->getComponent(0).size());
		cout << data.size() <<' ';
	}
	autocorr ac;
	autocorrInfo info;
	info = ac.autocorrelation(data, false);
	cout << "With Q = " << q << ":" << endl << info << endl;
	system("PAUSE");
	
	data.clear(); //reset data vector

	//
	knot->getComponent(0).setZ(z);
	initialComp0.readFromCoords(is);
	//knot = new clkConformationBfacf3(initialComp0);
	knot->step(1000000);
	for (int i = 0; i < 10000; i++){
		knot->step(100000);
		data.push_back(knot->getComponent(0).size());
	}
	info = ac.autocorrelation(data, false);
	cout << "Without Q:" << endl << info;
	

	/*
	//writhe sampling test code: currently testing the unknot
clkConformationBfacf3* knot; 
	vector<double> data;
	ifstream is;
	ofstream out;
	out.open("results/writhe_histogram_trefoil");
	is.open("initial/3_1");
	clkConformationAsList initialComp0;
	initialComp0.readFromCoords(is);
	knot = new clkConformationBfacf3(initialComp0);
	//knot->getComponent(0).setZ(.19273);
	double z = .2000;
	int q = 3;
	knot->setSeed(137); //set seed for test purposes
	//knot->init_Q(knot->getComponent(0).size(), knot->getZ(),3,false);
	knot->stepQ(1000000000,q,z);
	for (int i = 0; i < 20000; i++){
		knot->stepQ(650000,q,z);
		double writhe, acn;
		clkConformationAsList convert(knot->getComponent(0));
		convert.writheACN(writhe, acn);
		out << writhe << endl;
		data.push_back(writhe);
		cout << data.size() <<endl;
	}
	out.close();
	autocorr ac;
	autocorrInfo info;
	info = ac.autocorrelation(data, false);
	cout << "With Q = " << q << ":" << endl << info << endl;

	data.clear(); //reset data vector

	//
	knot->getComponent(0).setZ(z);
	initialComp0.readFromCoords(is);
	//knot = new clkConformationBfacf3(initialComp0);
	knot->step(1000000);
	for (int i = 0; i < 5000; i++){
		knot->step(25000);
		data.push_back(knot->getComponent(0).size());
	}
	info = ac.autocorrelation(data, false);
	cout << "Without Q:" << endl << info;
	*/
   	_CrtDumpMemoryLeaks();
	system("PAUSE");
	return 0;
}; 

