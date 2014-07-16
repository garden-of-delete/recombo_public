#include <testFramework.h>

#include <vector>
#include <string>

#include <threevector.h>
#include <genericConformation.h>
#include <clkCigar.h>
#include <clkConformationAsList.h>
#include <clkConformationBfacf3.h>
#include <bfacfProbabilities.h>
#include <bfacfProbabilitiesFromZ.h>
#include <bfacfProbabilitiesFromZFixed.h>
#include <pseudorandom.h>

#include <sstream>

#include <algorithm>
#include <deque>
#include <list>
#include <cmath>

#define NUMBER_OF_TREFOILS 5
#define TREFOIL_LENGTH 100
#define EPSILON 1.0e-10
#define DELTA 1.0e-4

using namespace std;

class testClkConformationBfacf3 : public clkConformationBfacf3{
public:
	testClkConformationBfacf3(const clkConformationAsList knot) : clkConformationBfacf3(knot) {};
	friend bool testClkConformationBfacf3::testPrecomputedBfacf3Probs();
};

class clktestclass
{
public:
	clktestclass();
	virtual ~clktestclass();
	void setUp();
	
	vector<double> expectedWr;
	vector<double> expectedAcn;
	vector<double> expectedRog;

	string manyTrefoilString;
	string manyTrefoilCoordString;

	vector<vector<int> > trefoilCoords;
	vector<vector<threevector<int> > > trefoilCoordsAsThreevectors;

	clk* clkTrefoil0;
	clk* clkTrefoil1;
	clk* clkTrefoil2;
	clk* clkTrefoil3;
	clk* clkTrefoil4;

};

clktestclass::clktestclass() : expectedAcn(NUMBER_OF_TREFOILS),
expectedRog(NUMBER_OF_TREFOILS), expectedWr(NUMBER_OF_TREFOILS),
trefoilCoords(NUMBER_OF_TREFOILS, vector<int>(3 * TREFOIL_LENGTH)),
trefoilCoordsAsThreevectors(NUMBER_OF_TREFOILS, vector<threevector<int> >(TREFOIL_LENGTH)) 
{
	setUp(); 
}

clktestclass::~clktestclass()
{
   delete clkTrefoil0;
   delete clkTrefoil1;
   delete clkTrefoil2;
   delete clkTrefoil3;
   delete clkTrefoil4;
}

void clktestclass::setUp()
{
   double arrayExpectedROG[] = {3.345594,
      3.949506,
      3.158006,
      3.358124,
      4.206043};
   double arrayExpectedWr[] = {3.750000,
      3.500000,
      1.750000,
      1.500000,
      3.500000};
   double arrayExpectedACN[] = {26.500493,
      18.401150,
      26.100593,
      25.527367,
      20.577284};
   copy(&arrayExpectedROG[0], &arrayExpectedROG[NUMBER_OF_TREFOILS], expectedRog.begin());
   copy(&arrayExpectedWr[0], &arrayExpectedWr[NUMBER_OF_TREFOILS], expectedWr.begin());
   copy(&arrayExpectedACN[0], &arrayExpectedACN[NUMBER_OF_TREFOILS], expectedAcn.begin());

   manyTrefoilString = " 20 5 -14 20 5 -15 21 5 -15 21 5 -16 21 4 -16 21 3 -16 22 3 -16 22 4 -16 22 4 -15 23 4 -15 24 4 -15 25 4 -15 25 4 -16 25 5 -16 25 5 -15 25 6 -15 24 6 -15 23 6 -15 23 5 -15 23 5 -14 22 5 -14 22 6 -14 22 6 -13 22 5 -13 21 5 -13 21 5 -12 20 5 -12 20 6 -12 20 6 -13 20 5 -13 20 4 -13 20 4 -14 20 3 -14 20 2 -14 19 2 -14 19 3 -14 19 4 -14 18 4 -14 18 4 -15 18 3 -15 19 3 -15 19 4 -15 19 5 -15 18 5 -15 18 5 -16 17 5 -16 17 6 -16 18 6 -16 19 6 -16 20 6 -16 20 5 -16 20 5 -17 20 5 -18 20 5 -19 20 6 -19 21 6 -19 22 6 -19 22 5 -19 23 5 -19 24 5 -19 24 5 -18 24 4 -18 24 4 -17 24 4 -16 24 5 -16 24 5 -15 24 5 -14 24 4 -14 24 4 -13 24 5 -13 24 6 -13 23 6 -13 23 6 -12 23 7 -12 23 8 -12 23 9 -12 23 9 -11 22 9 -11 21 9 -11 21 8 -11 21 8 -10 20 8 -10 20 8 -11 20 8 -12 20 8 -13 19 8 -13 19 8 -12 19 7 -12 19 7 -13 20 7 -13 21 7 -13 21 7 -14 21 7 -15 20 7 -15 20 8 -15 19 8 -15 19 7 -15 19 6 -15 19 6 -14 19 5 -14\n\n 20 -7 -23 19 -7 -23 19 -7 -22 19 -7 -21 19 -6 -21 19 -6 -20 18 -6 -20 18 -6 -19 18 -5 -19 19 -5 -19 19 -5 -18 19 -5 -17 19 -5 -16 20 -5 -16 20 -5 -17 21 -5 -17 21 -4 -17 21 -3 -17 22 -3 -17 22 -3 -18 22 -2 -18 22 -2 -19 23 -2 -19 23 -2 -20 24 -2 -20 25 -2 -20 25 -3 -20 25 -4 -20 26 -4 -20 26 -5 -20 26 -5 -19 25 -5 -19 25 -6 -19 25 -7 -19 25 -7 -18 24 -7 -18 24 -7 -19 24 -8 -19 24 -9 -19 23 -9 -19 23 -9 -20 23 -9 -21 22 -9 -21 22 -8 -21 22 -8 -20 21 -8 -20 20 -8 -20 20 -8 -21 20 -7 -21 20 -6 -21 20 -6 -22 20 -5 -22 19 -5 -22 18 -5 -22 17 -5 -22 17 -6 -22 17 -6 -21 17 -7 -21 18 -7 -21 18 -7 -20 18 -7 -19 17 -7 -19 17 -7 -18 17 -8 -18 17 -8 -17 17 -8 -16 17 -9 -16 17 -9 -15 18 -9 -15 18 -8 -15 18 -7 -15 19 -7 -15 19 -7 -14 20 -7 -14 21 -7 -14 21 -8 -14 21 -8 -15 22 -8 -15 22 -8 -16 22 -8 -17 22 -7 -17 22 -7 -18 22 -6 -18 22 -5 -18 21 -5 -18 21 -5 -19 20 -5 -19 20 -5 -20 20 -6 -20 21 -6 -20 22 -6 -20 22 -6 -21 22 -6 -22 22 -6 -23 22 -6 -24 21 -6 -24 21 -6 -25 20 -6 -25 20 -6 -24 20 -6 -23\n\n 14 -5 -17 15 -5 -17 16 -5 -17 16 -6 -17 16 -7 -17 16 -8 -17 17 -8 -17 17 -8 -18 16 -8 -18 15 -8 -18 15 -8 -19 15 -8 -20 14 -8 -20 14 -7 -20 13 -7 -20 13 -6 -20 12 -6 -20 12 -6 -19 12 -6 -18 12 -6 -17 12 -6 -16 11 -6 -16 11 -6 -15 11 -7 -15 11 -7 -16 10 -7 -16 10 -8 -16 11 -8 -16 11 -8 -17 12 -8 -17 12 -7 -17 12 -7 -16 13 -7 -16 13 -7 -15 13 -6 -15 13 -6 -14 13 -6 -13 14 -6 -13 14 -6 -14 14 -6 -15 14 -7 -15 14 -7 -14 14 -8 -14 14 -8 -15 14 -8 -16 14 -9 -16 14 -9 -17 15 -9 -17 15 -8 -17 15 -7 -17 15 -7 -18 15 -7 -19 15 -6 -19 15 -5 -19 15 -5 -18 14 -5 -18 14 -4 -18 14 -3 -18 14 -3 -17 14 -3 -16 14 -2 -16 13 -2 -16 12 -2 -16 12 -3 -16 12 -3 -17 12 -3 -18 11 -3 -18 11 -4 -18 10 -4 -18 10 -4 -19 11 -4 -19 12 -4 -19 12 -3 -19 12 -3 -20 13 -3 -20 14 -3 -20 14 -4 -20 14 -4 -21 13 -4 -21 13 -5 -21 13 -6 -21 12 -6 -21 12 -7 -21 11 -7 -21 11 -7 -20 11 -8 -20 12 -8 -20 12 -8 -19 12 -7 -19 13 -7 -19 14 -7 -19 14 -8 -19 14 -9 -19 15 -9 -19 15 -9 -18 14 -9 -18 14 -8 -18 14 -8 -17 14 -7 -17 14 -6 -17\n\n 13 -1 -16 12 -1 -16 12 -2 -16 12 -3 -16 12 -4 -16 11 -4 -16 11 -4 -15 11 -5 -15 12 -5 -15 12 -5 -14 11 -5 -14 11 -6 -14 11 -7 -14 11 -7 -15 11 -7 -16 11 -7 -17 11 -8 -17 11 -9 -17 12 -9 -17 12 -8 -17 12 -7 -17 13 -7 -17 13 -7 -18 13 -6 -18 13 -6 -19 12 -6 -19 11 -6 -19 11 -6 -18 11 -6 -17 11 -5 -17 11 -5 -18 12 -5 -18 13 -5 -18 14 -5 -18 14 -6 -18 14 -7 -18 14 -8 -18 14 -8 -19 14 -7 -19 14 -6 -19 14 -5 -19 14 -5 -20 15 -5 -20 15 -5 -19 15 -4 -19 14 -4 -19 13 -4 -19 13 -4 -18 12 -4 -18 11 -4 -18 10 -4 -18 10 -3 -18 9 -3 -18 8 -3 -18 8 -3 -19 8 -4 -19 9 -4 -19 9 -3 -19 10 -3 -19 10 -3 -20 10 -2 -20 10 -2 -21 11 -2 -21 12 -2 -21 13 -2 -21 13 -1 -21 12 -1 -21 12 -1 -20 11 -1 -20 11 0 -20 11 1 -20 12 1 -20 12 0 -20 13 0 -20 13 -1 -20 13 -1 -19 14 -1 -19 14 -1 -20 15 -1 -20 15 -2 -20 14 -2 -20 14 -2 -19 13 -2 -19 13 -2 -18 13 -2 -17 13 -2 -16 14 -2 -16 14 -3 -16 13 -3 -16 13 -3 -15 12 -3 -15 11 -3 -15 11 -3 -16 11 -3 -17 12 -3 -17 13 -3 -17 14 -3 -17 14 -2 -17 14 -1 -17 14 -1 -16\n\n 19 5 -17 20 5 -17 20 6 -17 19 6 -17 19 6 -16 18 6 -16 18 6 -15 18 5 -15 18 5 -14 18 6 -14 18 7 -14 18 7 -13 18 7 -12 18 7 -11 18 6 -11 18 6 -12 17 6 -12 17 5 -12 17 4 -12 17 3 -12 17 3 -13 17 2 -13 17 2 -14 17 2 -15 16 2 -15 16 2 -16 16 3 -16 16 4 -16 16 5 -16 17 5 -16 18 5 -16 19 5 -16 19 5 -15 20 5 -15 20 5 -14 21 5 -14 21 6 -14 21 6 -15 21 6 -16 22 6 -16 22 6 -17 22 7 -17 22 7 -16 21 7 -16 20 7 -16 20 7 -15 20 7 -14 20 7 -13 19 7 -13 19 8 -13 18 8 -13 18 8 -12 18 9 -12 18 9 -13 17 9 -13 17 9 -12 16 9 -12 16 10 -12 15 10 -12 15 10 -13 14 10 -13 14 9 -13 14 8 -13 14 8 -14 13 8 -14 13 9 -14 13 10 -14 12 10 -14 11 10 -14 11 9 -14 11 9 -15 10 9 -15 10 9 -16 10 9 -17 10 8 -17 10 7 -17 11 7 -17 12 7 -17 12 8 -17 12 8 -16 11 8 -16 11 7 -16 11 7 -15 11 6 -15 12 6 -15 13 6 -15 14 6 -15 14 5 -15 15 5 -15 15 5 -14 16 5 -14 16 5 -15 16 4 -15 16 4 -14 17 4 -14 17 4 -15 18 4 -15 19 4 -15 19 4 -16 19 4 -17";   
  
   double trefoil0[] = {20, 5, -14, 20, 5, -15, 21, 5, -15, 21, 5, -16, 21, 4, -16, 21, 3, -16, 22, 3, -16, 22, 4, -16, 22, 4, -15, 23, 4, -15, 24, 4, -15, 25, 4, -15, 25, 4, -16, 25, 5, -16, 25, 5, -15, 25, 6, -15, 24, 6, -15, 23, 6, -15, 23, 5, -15, 23, 5, -14, 22, 5, -14, 22, 6, -14, 22, 6, -13, 22, 5, -13, 21, 5, -13, 21, 5, -12, 20, 5, -12, 20, 6, -12, 20, 6, -13, 20, 5, -13, 20, 4, -13, 20, 4, -14, 20, 3, -14, 20, 2, -14, 19, 2, -14, 19, 3, -14, 19, 4, -14, 18, 4, -14, 18, 4, -15, 18, 3, -15, 19, 3, -15, 19, 4, -15, 19, 5, -15, 18, 5, -15, 18, 5, -16, 17, 5, -16, 17, 6, -16, 18, 6, -16, 19, 6, -16, 20, 6, -16, 20, 5, -16, 20, 5, -17, 20, 5, -18, 20, 5, -19, 20, 6, -19, 21, 6, -19, 22, 6, -19, 22, 5, -19, 23, 5, -19, 24, 5, -19, 24, 5, -18, 24, 4, -18, 24, 4, -17, 24, 4, -16, 24, 5, -16, 24, 5, -15, 24, 5, -14, 24, 4, -14, 24, 4, -13, 24, 5, -13, 24, 6, -13, 23, 6, -13, 23, 6, -12, 23, 7, -12, 23, 8, -12, 23, 9, -12, 23, 9, -11, 22, 9, -11, 21, 9, -11, 21, 8, -11, 21, 8, -10, 20, 8, -10, 20, 8, -11, 20, 8, -12, 20, 8, -13, 19, 8, -13, 19, 8, -12, 19, 7, -12, 19, 7, -13, 20, 7, -13, 21, 7, -13, 21, 7, -14, 21, 7, -15, 20, 7, -15, 20, 8, -15, 19, 8, -15, 19, 7, -15, 19, 6, -15, 19, 6, -14, 19, 5, -14};
   double trefoil1[] = {20, -7, -23, 19, -7, -23, 19, -7, -22, 19, -7, -21, 19, -6, -21, 19, -6, -20, 18, -6, -20, 18, -6, -19, 18, -5, -19, 19, -5, -19, 19, -5, -18, 19, -5, -17, 19, -5, -16, 20, -5, -16, 20, -5, -17, 21, -5, -17, 21, -4, -17, 21, -3, -17, 22, -3, -17, 22, -3, -18, 22, -2, -18, 22, -2, -19, 23, -2, -19, 23, -2, -20, 24, -2, -20, 25, -2, -20, 25, -3, -20, 25, -4, -20, 26, -4, -20, 26, -5, -20, 26, -5, -19, 25, -5, -19, 25, -6, -19, 25, -7, -19, 25, -7, -18, 24, -7, -18, 24, -7, -19, 24, -8, -19, 24, -9, -19, 23, -9, -19, 23, -9, -20, 23, -9, -21, 22, -9, -21, 22, -8, -21, 22, -8, -20, 21, -8, -20, 20, -8, -20, 20, -8, -21, 20, -7, -21, 20, -6, -21, 20, -6, -22, 20, -5, -22, 19, -5, -22, 18, -5, -22, 17, -5, -22, 17, -6, -22, 17, -6, -21, 17, -7, -21, 18, -7, -21, 18, -7, -20, 18, -7, -19, 17, -7, -19, 17, -7, -18, 17, -8, -18, 17, -8, -17, 17, -8, -16, 17, -9, -16, 17, -9, -15, 18, -9, -15, 18, -8, -15, 18, -7, -15, 19, -7, -15, 19, -7, -14, 20, -7, -14, 21, -7, -14, 21, -8, -14, 21, -8, -15, 22, -8, -15, 22, -8, -16, 22, -8, -17, 22, -7, -17, 22, -7, -18, 22, -6, -18, 22, -5, -18, 21, -5, -18, 21, -5, -19, 20, -5, -19, 20, -5, -20, 20, -6, -20, 21, -6, -20, 22, -6, -20, 22, -6, -21, 22, -6, -22, 22, -6, -23, 22, -6, -24, 21, -6, -24, 21, -6, -25, 20, -6, -25, 20, -6, -24, 20, -6, -23};
   double trefoil2[] = {14, -5, -17, 15, -5, -17, 16, -5, -17, 16, -6, -17, 16, -7, -17, 16, -8, -17, 17, -8, -17, 17, -8, -18, 16, -8, -18, 15, -8, -18, 15, -8, -19, 15, -8, -20, 14, -8, -20, 14, -7, -20, 13, -7, -20, 13, -6, -20, 12, -6, -20, 12, -6, -19, 12, -6, -18, 12, -6, -17, 12, -6, -16, 11, -6, -16, 11, -6, -15, 11, -7, -15, 11, -7, -16, 10, -7, -16, 10, -8, -16, 11, -8, -16, 11, -8, -17, 12, -8, -17, 12, -7, -17, 12, -7, -16, 13, -7, -16, 13, -7, -15, 13, -6, -15, 13, -6, -14, 13, -6, -13, 14, -6, -13, 14, -6, -14, 14, -6, -15, 14, -7, -15, 14, -7, -14, 14, -8, -14, 14, -8, -15, 14, -8, -16, 14, -9, -16, 14, -9, -17, 15, -9, -17, 15, -8, -17, 15, -7, -17, 15, -7, -18, 15, -7, -19, 15, -6, -19, 15, -5, -19, 15, -5, -18, 14, -5, -18, 14, -4, -18, 14, -3, -18, 14, -3, -17, 14, -3, -16, 14, -2, -16, 13, -2, -16, 12, -2, -16, 12, -3, -16, 12, -3, -17, 12, -3, -18, 11, -3, -18, 11, -4, -18, 10, -4, -18, 10, -4, -19, 11, -4, -19, 12, -4, -19, 12, -3, -19, 12, -3, -20, 13, -3, -20, 14, -3, -20, 14, -4, -20, 14, -4, -21, 13, -4, -21, 13, -5, -21, 13, -6, -21, 12, -6, -21, 12, -7, -21, 11, -7, -21, 11, -7, -20, 11, -8, -20, 12, -8, -20, 12, -8, -19, 12, -7, -19, 13, -7, -19, 14, -7, -19, 14, -8, -19, 14, -9, -19, 15, -9, -19, 15, -9, -18, 14, -9, -18, 14, -8, -18, 14, -8, -17, 14, -7, -17, 14, -6, -17};
   double trefoil3[] = {13, -1, -16, 12, -1, -16, 12, -2, -16, 12, -3, -16, 12, -4, -16, 11, -4, -16, 11, -4, -15, 11, -5, -15, 12, -5, -15, 12, -5, -14, 11, -5, -14, 11, -6, -14, 11, -7, -14, 11, -7, -15, 11, -7, -16, 11, -7, -17, 11, -8, -17, 11, -9, -17, 12, -9, -17, 12, -8, -17, 12, -7, -17, 13, -7, -17, 13, -7, -18, 13, -6, -18, 13, -6, -19, 12, -6, -19, 11, -6, -19, 11, -6, -18, 11, -6, -17, 11, -5, -17, 11, -5, -18, 12, -5, -18, 13, -5, -18, 14, -5, -18, 14, -6, -18, 14, -7, -18, 14, -8, -18, 14, -8, -19, 14, -7, -19, 14, -6, -19, 14, -5, -19, 14, -5, -20, 15, -5, -20, 15, -5, -19, 15, -4, -19, 14, -4, -19, 13, -4, -19, 13, -4, -18, 12, -4, -18, 11, -4, -18, 10, -4, -18, 10, -3, -18, 9, -3, -18, 8, -3, -18, 8, -3, -19, 8, -4, -19, 9, -4, -19, 9, -3, -19, 10, -3, -19, 10, -3, -20, 10, -2, -20, 10, -2, -21, 11, -2, -21, 12, -2, -21, 13, -2, -21, 13, -1, -21, 12, -1, -21, 12, -1, -20, 11, -1, -20, 11, 0, -20, 11, 1, -20, 12, 1, -20, 12, 0, -20, 13, 0, -20, 13, -1, -20, 13, -1, -19, 14, -1, -19, 14, -1, -20, 15, -1, -20, 15, -2, -20, 14, -2, -20, 14, -2, -19, 13, -2, -19, 13, -2, -18, 13, -2, -17, 13, -2, -16, 14, -2, -16, 14, -3, -16, 13, -3, -16, 13, -3, -15, 12, -3, -15, 11, -3, -15, 11, -3, -16, 11, -3, -17, 12, -3, -17, 13, -3, -17, 14, -3, -17, 14, -2, -17, 14, -1, -17, 14, -1, -16};
   double trefoil4[] = {19, 5, -17, 20, 5, -17, 20, 6, -17, 19, 6, -17, 19, 6, -16, 18, 6, -16, 18, 6, -15, 18, 5, -15, 18, 5, -14, 18, 6, -14, 18, 7, -14, 18, 7, -13, 18, 7, -12, 18, 7, -11, 18, 6, -11, 18, 6, -12, 17, 6, -12, 17, 5, -12, 17, 4, -12, 17, 3, -12, 17, 3, -13, 17, 2, -13, 17, 2, -14, 17, 2, -15, 16, 2, -15, 16, 2, -16, 16, 3, -16, 16, 4, -16, 16, 5, -16, 17, 5, -16, 18, 5, -16, 19, 5, -16, 19, 5, -15, 20, 5, -15, 20, 5, -14, 21, 5, -14, 21, 6, -14, 21, 6, -15, 21, 6, -16, 22, 6, -16, 22, 6, -17, 22, 7, -17, 22, 7, -16, 21, 7, -16, 20, 7, -16, 20, 7, -15, 20, 7, -14, 20, 7, -13, 19, 7, -13, 19, 8, -13, 18, 8, -13, 18, 8, -12, 18, 9, -12, 18, 9, -13, 17, 9, -13, 17, 9, -12, 16, 9, -12, 16, 10, -12, 15, 10, -12, 15, 10, -13, 14, 10, -13, 14, 9, -13, 14, 8, -13, 14, 8, -14, 13, 8, -14, 13, 9, -14, 13, 10, -14, 12, 10, -14, 11, 10, -14, 11, 9, -14, 11, 9, -15, 10, 9, -15, 10, 9, -16, 10, 9, -17, 10, 8, -17, 10, 7, -17, 11, 7, -17, 12, 7, -17, 12, 8, -17, 12, 8, -16, 11, 8, -16, 11, 7, -16, 11, 7, -15, 11, 6, -15, 12, 6, -15, 13, 6, -15, 14, 6, -15, 14, 5, -15, 15, 5, -15, 15, 5, -14, 16, 5, -14, 16, 5, -15, 16, 4, -15, 16, 4, -14, 17, 4, -14, 17, 4, -15, 18, 4, -15, 19, 4, -15, 19, 4, -16, 19, 4, -17};

   copy(&trefoil0[0], &trefoil0[3 * TREFOIL_LENGTH], trefoilCoords[0].begin());
   copy(&trefoil1[0], &trefoil1[3 * TREFOIL_LENGTH], trefoilCoords[1].begin());
   copy(&trefoil2[0], &trefoil2[3 * TREFOIL_LENGTH], trefoilCoords[2].begin());
   copy(&trefoil3[0], &trefoil3[3 * TREFOIL_LENGTH], trefoilCoords[3].begin());
   copy(&trefoil4[0], &trefoil4[3 * TREFOIL_LENGTH], trefoilCoords[4].begin());

   for (int i = 0; i < TREFOIL_LENGTH; i++)
   {
      trefoilCoordsAsThreevectors[0][i].set(trefoil0[3 * i + 0], trefoil0[3 * i + 1], trefoil0[3 * i + 2]);
      trefoilCoordsAsThreevectors[1][i].set(trefoil1[3 * i + 0], trefoil1[3 * i + 1], trefoil1[3 * i + 2]);
      trefoilCoordsAsThreevectors[2][i].set(trefoil2[3 * i + 0], trefoil2[3 * i + 1], trefoil2[3 * i + 2]);
      trefoilCoordsAsThreevectors[3][i].set(trefoil3[3 * i + 0], trefoil3[3 * i + 1], trefoil3[3 * i + 2]);
      trefoilCoordsAsThreevectors[4][i].set(trefoil4[3 * i + 0], trefoil4[3 * i + 1], trefoil4[3 * i + 2]);
   }

   int itrefoil0[] = {20, 5, -14, 20, 5, -15, 21, 5, -15, 21, 5, -16, 21, 4, -16, 21, 3, -16, 22, 3, -16, 22, 4, -16, 22, 4, -15, 23, 4, -15, 24, 4, -15, 25, 4, -15, 25, 4, -16, 25, 5, -16, 25, 5, -15, 25, 6, -15, 24, 6, -15, 23, 6, -15, 23, 5, -15, 23, 5, -14, 22, 5, -14, 22, 6, -14, 22, 6, -13, 22, 5, -13, 21, 5, -13, 21, 5, -12, 20, 5, -12, 20, 6, -12, 20, 6, -13, 20, 5, -13, 20, 4, -13, 20, 4, -14, 20, 3, -14, 20, 2, -14, 19, 2, -14, 19, 3, -14, 19, 4, -14, 18, 4, -14, 18, 4, -15, 18, 3, -15, 19, 3, -15, 19, 4, -15, 19, 5, -15, 18, 5, -15, 18, 5, -16, 17, 5, -16, 17, 6, -16, 18, 6, -16, 19, 6, -16, 20, 6, -16, 20, 5, -16, 20, 5, -17, 20, 5, -18, 20, 5, -19, 20, 6, -19, 21, 6, -19, 22, 6, -19, 22, 5, -19, 23, 5, -19, 24, 5, -19, 24, 5, -18, 24, 4, -18, 24, 4, -17, 24, 4, -16, 24, 5, -16, 24, 5, -15, 24, 5, -14, 24, 4, -14, 24, 4, -13, 24, 5, -13, 24, 6, -13, 23, 6, -13, 23, 6, -12, 23, 7, -12, 23, 8, -12, 23, 9, -12, 23, 9, -11, 22, 9, -11, 21, 9, -11, 21, 8, -11, 21, 8, -10, 20, 8, -10, 20, 8, -11, 20, 8, -12, 20, 8, -13, 19, 8, -13, 19, 8, -12, 19, 7, -12, 19, 7, -13, 20, 7, -13, 21, 7, -13, 21, 7, -14, 21, 7, -15, 20, 7, -15, 20, 8, -15, 19, 8, -15, 19, 7, -15, 19, 6, -15, 19, 6, -14, 19, 5, -14};
   int itrefoil1[] = {20, -7, -23, 19, -7, -23, 19, -7, -22, 19, -7, -21, 19, -6, -21, 19, -6, -20, 18, -6, -20, 18, -6, -19, 18, -5, -19, 19, -5, -19, 19, -5, -18, 19, -5, -17, 19, -5, -16, 20, -5, -16, 20, -5, -17, 21, -5, -17, 21, -4, -17, 21, -3, -17, 22, -3, -17, 22, -3, -18, 22, -2, -18, 22, -2, -19, 23, -2, -19, 23, -2, -20, 24, -2, -20, 25, -2, -20, 25, -3, -20, 25, -4, -20, 26, -4, -20, 26, -5, -20, 26, -5, -19, 25, -5, -19, 25, -6, -19, 25, -7, -19, 25, -7, -18, 24, -7, -18, 24, -7, -19, 24, -8, -19, 24, -9, -19, 23, -9, -19, 23, -9, -20, 23, -9, -21, 22, -9, -21, 22, -8, -21, 22, -8, -20, 21, -8, -20, 20, -8, -20, 20, -8, -21, 20, -7, -21, 20, -6, -21, 20, -6, -22, 20, -5, -22, 19, -5, -22, 18, -5, -22, 17, -5, -22, 17, -6, -22, 17, -6, -21, 17, -7, -21, 18, -7, -21, 18, -7, -20, 18, -7, -19, 17, -7, -19, 17, -7, -18, 17, -8, -18, 17, -8, -17, 17, -8, -16, 17, -9, -16, 17, -9, -15, 18, -9, -15, 18, -8, -15, 18, -7, -15, 19, -7, -15, 19, -7, -14, 20, -7, -14, 21, -7, -14, 21, -8, -14, 21, -8, -15, 22, -8, -15, 22, -8, -16, 22, -8, -17, 22, -7, -17, 22, -7, -18, 22, -6, -18, 22, -5, -18, 21, -5, -18, 21, -5, -19, 20, -5, -19, 20, -5, -20, 20, -6, -20, 21, -6, -20, 22, -6, -20, 22, -6, -21, 22, -6, -22, 22, -6, -23, 22, -6, -24, 21, -6, -24, 21, -6, -25, 20, -6, -25, 20, -6, -24, 20, -6, -23};
   int itrefoil2[] = {14, -5, -17, 15, -5, -17, 16, -5, -17, 16, -6, -17, 16, -7, -17, 16, -8, -17, 17, -8, -17, 17, -8, -18, 16, -8, -18, 15, -8, -18, 15, -8, -19, 15, -8, -20, 14, -8, -20, 14, -7, -20, 13, -7, -20, 13, -6, -20, 12, -6, -20, 12, -6, -19, 12, -6, -18, 12, -6, -17, 12, -6, -16, 11, -6, -16, 11, -6, -15, 11, -7, -15, 11, -7, -16, 10, -7, -16, 10, -8, -16, 11, -8, -16, 11, -8, -17, 12, -8, -17, 12, -7, -17, 12, -7, -16, 13, -7, -16, 13, -7, -15, 13, -6, -15, 13, -6, -14, 13, -6, -13, 14, -6, -13, 14, -6, -14, 14, -6, -15, 14, -7, -15, 14, -7, -14, 14, -8, -14, 14, -8, -15, 14, -8, -16, 14, -9, -16, 14, -9, -17, 15, -9, -17, 15, -8, -17, 15, -7, -17, 15, -7, -18, 15, -7, -19, 15, -6, -19, 15, -5, -19, 15, -5, -18, 14, -5, -18, 14, -4, -18, 14, -3, -18, 14, -3, -17, 14, -3, -16, 14, -2, -16, 13, -2, -16, 12, -2, -16, 12, -3, -16, 12, -3, -17, 12, -3, -18, 11, -3, -18, 11, -4, -18, 10, -4, -18, 10, -4, -19, 11, -4, -19, 12, -4, -19, 12, -3, -19, 12, -3, -20, 13, -3, -20, 14, -3, -20, 14, -4, -20, 14, -4, -21, 13, -4, -21, 13, -5, -21, 13, -6, -21, 12, -6, -21, 12, -7, -21, 11, -7, -21, 11, -7, -20, 11, -8, -20, 12, -8, -20, 12, -8, -19, 12, -7, -19, 13, -7, -19, 14, -7, -19, 14, -8, -19, 14, -9, -19, 15, -9, -19, 15, -9, -18, 14, -9, -18, 14, -8, -18, 14, -8, -17, 14, -7, -17, 14, -6, -17};
   int itrefoil3[] = {13, -1, -16, 12, -1, -16, 12, -2, -16, 12, -3, -16, 12, -4, -16, 11, -4, -16, 11, -4, -15, 11, -5, -15, 12, -5, -15, 12, -5, -14, 11, -5, -14, 11, -6, -14, 11, -7, -14, 11, -7, -15, 11, -7, -16, 11, -7, -17, 11, -8, -17, 11, -9, -17, 12, -9, -17, 12, -8, -17, 12, -7, -17, 13, -7, -17, 13, -7, -18, 13, -6, -18, 13, -6, -19, 12, -6, -19, 11, -6, -19, 11, -6, -18, 11, -6, -17, 11, -5, -17, 11, -5, -18, 12, -5, -18, 13, -5, -18, 14, -5, -18, 14, -6, -18, 14, -7, -18, 14, -8, -18, 14, -8, -19, 14, -7, -19, 14, -6, -19, 14, -5, -19, 14, -5, -20, 15, -5, -20, 15, -5, -19, 15, -4, -19, 14, -4, -19, 13, -4, -19, 13, -4, -18, 12, -4, -18, 11, -4, -18, 10, -4, -18, 10, -3, -18, 9, -3, -18, 8, -3, -18, 8, -3, -19, 8, -4, -19, 9, -4, -19, 9, -3, -19, 10, -3, -19, 10, -3, -20, 10, -2, -20, 10, -2, -21, 11, -2, -21, 12, -2, -21, 13, -2, -21, 13, -1, -21, 12, -1, -21, 12, -1, -20, 11, -1, -20, 11, 0, -20, 11, 1, -20, 12, 1, -20, 12, 0, -20, 13, 0, -20, 13, -1, -20, 13, -1, -19, 14, -1, -19, 14, -1, -20, 15, -1, -20, 15, -2, -20, 14, -2, -20, 14, -2, -19, 13, -2, -19, 13, -2, -18, 13, -2, -17, 13, -2, -16, 14, -2, -16, 14, -3, -16, 13, -3, -16, 13, -3, -15, 12, -3, -15, 11, -3, -15, 11, -3, -16, 11, -3, -17, 12, -3, -17, 13, -3, -17, 14, -3, -17, 14, -2, -17, 14, -1, -17, 14, -1, -16};
   int itrefoil4[] = {19, 5, -17, 20, 5, -17, 20, 6, -17, 19, 6, -17, 19, 6, -16, 18, 6, -16, 18, 6, -15, 18, 5, -15, 18, 5, -14, 18, 6, -14, 18, 7, -14, 18, 7, -13, 18, 7, -12, 18, 7, -11, 18, 6, -11, 18, 6, -12, 17, 6, -12, 17, 5, -12, 17, 4, -12, 17, 3, -12, 17, 3, -13, 17, 2, -13, 17, 2, -14, 17, 2, -15, 16, 2, -15, 16, 2, -16, 16, 3, -16, 16, 4, -16, 16, 5, -16, 17, 5, -16, 18, 5, -16, 19, 5, -16, 19, 5, -15, 20, 5, -15, 20, 5, -14, 21, 5, -14, 21, 6, -14, 21, 6, -15, 21, 6, -16, 22, 6, -16, 22, 6, -17, 22, 7, -17, 22, 7, -16, 21, 7, -16, 20, 7, -16, 20, 7, -15, 20, 7, -14, 20, 7, -13, 19, 7, -13, 19, 8, -13, 18, 8, -13, 18, 8, -12, 18, 9, -12, 18, 9, -13, 17, 9, -13, 17, 9, -12, 16, 9, -12, 16, 10, -12, 15, 10, -12, 15, 10, -13, 14, 10, -13, 14, 9, -13, 14, 8, -13, 14, 8, -14, 13, 8, -14, 13, 9, -14, 13, 10, -14, 12, 10, -14, 11, 10, -14, 11, 9, -14, 11, 9, -15, 10, 9, -15, 10, 9, -16, 10, 9, -17, 10, 8, -17, 10, 7, -17, 11, 7, -17, 12, 7, -17, 12, 8, -17, 12, 8, -16, 11, 8, -16, 11, 7, -16, 11, 7, -15, 11, 6, -15, 12, 6, -15, 13, 6, -15, 14, 6, -15, 14, 5, -15, 15, 5, -15, 15, 5, -14, 16, 5, -14, 16, 5, -15, 16, 4, -15, 16, 4, -14, 17, 4, -14, 17, 4, -15, 18, 4, -15, 19, 4, -15, 19, 4, -16, 19, 4, -17};

   clkTrefoil0 = new clkConformationAsList(itrefoil0, TREFOIL_LENGTH);
   clkTrefoil1 = new clkConformationAsList(itrefoil1, TREFOIL_LENGTH);
   clkTrefoil2 = new clkConformationAsList(itrefoil2, TREFOIL_LENGTH);
   clkTrefoil3 = new clkConformationAsList(itrefoil3, TREFOIL_LENGTH);
   clkTrefoil4 = new clkConformationAsList(itrefoil4, TREFOIL_LENGTH);
   
   stringstream s0, s1, s2, s3, s4;
   for(int i = 0; i < TREFOIL_LENGTH; i++)
   {
      s0 << itrefoil0[3*i] << " " << itrefoil0[3*i+1] << " " << itrefoil0[3*i+2] << endl;
      s1 << itrefoil1[3*i] << " " << itrefoil1[3*i+1] << " " << itrefoil1[3*i+2] << endl;
      s2 << itrefoil2[3*i] << " " << itrefoil2[3*i+1] << " " << itrefoil2[3*i+2] << endl;
      s3 << itrefoil3[3*i] << " " << itrefoil3[3*i+1] << " " << itrefoil3[3*i+2] << endl;
      s4 << itrefoil4[3*i] << " " << itrefoil4[3*i+1] << " " << itrefoil4[3*i+2] << endl;
   }
   s0 << "0" << endl;
   s1 << "1" << endl;
   s2 << "2" << endl;
   s3 << "3" << endl;
   s4 << "4" << endl;
   manyTrefoilCoordString += s0.str();
   manyTrefoilCoordString += s1.str();
   manyTrefoilCoordString += s2.str();
   manyTrefoilCoordString += s3.str();
   manyTrefoilCoordString += s4.str();
}

bool vectorsEqual(const vector<threevector<int> >& u, const vector<threevector<int> >& v)
{
   if (u.size() != v.size()) return false;

   for (int i = 0; i < u.size(); i++)
      if (u[i] != v[i]) return false;

   return true;
}

bool vectorsEqual(const vector<int>& u, const vector<threevector<int> >& v)
{
   if (u.size() != 3 * v.size()) return false;

   for (int i = 0; i < v.size(); i++)
      if (v[i] != threevector<int>(u[3 * i], u[3 * i + 1], u[3 * i + 2])) return false;

   return true;
}

bool testData()
{
	clktestclass clk;
	ASSERT_MESSAGE("Wrong length for expectedAcn.", clk.expectedAcn.size() == NUMBER_OF_TREFOILS);
	vector<double>::const_iterator i;

	ASSERT_MESSAGE("Wrong length for expectedWr.", clk.expectedRog.size() == NUMBER_OF_TREFOILS);
	ASSERT_MESSAGE("Wrong length for expectedRog.", clk.expectedWr.size() == NUMBER_OF_TREFOILS);

	ASSERT_MESSAGE("Wrong number of trefoils.", clk.trefoilCoords.size() == NUMBER_OF_TREFOILS);
	for (vector<vector<int> >::const_iterator j = clk.trefoilCoords.begin();
	   j != clk.trefoilCoords.end(); j++)
	{
		ASSERT_MESSAGE("Wrong length of trefoil.", j->size() == 3 * TREFOIL_LENGTH);
	}

	ASSERT_MESSAGE("Wrong number of trefoils as threevectors.", clk.trefoilCoordsAsThreevectors.size() == NUMBER_OF_TREFOILS);
	for (int j = 0; j < NUMBER_OF_TREFOILS; j++)
	{
		ASSERT_MESSAGE("Wrong length of trefoil as threevectors.", clk.trefoilCoordsAsThreevectors[j].size() == TREFOIL_LENGTH);

		ASSERT_MESSAGE("Trefoil as threevectors do not have expected coords.", vectorsEqual(clk.trefoilCoords[j], clk.trefoilCoordsAsThreevectors[j]));
	}

	return true;
}

bool testCopyScalarToVector(const vector<int>& u)
{
   // Copying a vector<int> to a vector<threevector<int> >
   int length = u.size() / 3;
   vector<threevector<int> > v(length);
   copyScalarToVector(u.begin(), u.end(), v.begin());

   ASSERT_MESSAGE("Copy does not match original.", vectorsEqual(u, v));

   // Copying a vector<int> to a deque<threevector<int> >.
   deque<threevector<int> > dq(length);
   copyScalarToVector(u.begin(), u.end(), dq.begin());

   ASSERT_MESSAGE("Deque copy does not match original in number of vertices size.", dq.size());
   ASSERT_MESSAGE("Deque copy does not match original.", equal(v.begin(), v.end(), dq.begin()));

   // Copying a vector<int> to a list<threevector<int> >.
   list<threevector<int> > l(length);
   copyScalarToVector(u.begin(), u.end(), l.begin());

   ASSERT_MESSAGE("List copy does not match original in number of vertices size.", l.size());
   ASSERT_MESSAGE("List copy does not match original.", equal(v.begin(), v.end(), l.begin()));

   return true;
}

bool testCopyScalarToVector()
{
	clktestclass clk;
 
   for (int i = 0; i < clk.trefoilCoords.size(); i++)
      if(!testCopyScalarToVector(clk.trefoilCoords[i])) return false;
   return true;
}

bool testGenericOutput(const vector<threevector<int> >& u)
{
   stringstream strforvector, strfordeque, strforlist;
   deque<threevector<int> > dq(u.size());
   list<threevector<int> > l(u.size());
   copy(u.begin(), u.end(), dq.begin());
   copy(u.begin(), u.end(), l.begin());

   outputVertices(u.begin(), u.end(), strforvector, "[", "]", "; ", "(", ")", ", ");
   outputVertices(dq.begin(), dq.end(), strfordeque, "[", "]", "; ", "(", ")", ", ");
   outputVertices(l.begin(), l.end(), strforlist, "[", "]", "; ", "(", ")", ", ");

   string expected(strforvector.str());
   ASSERT(expected == strfordeque.str());
   ASSERT(expected == strforlist.str());
   
   return true;
}

bool testGenericOutput()
{
	clktestclass clk;
 
   stringstream str1, str2, str3, str4;
   string expected1("1 2 3"), expected2("(1, 2, 3)"),
           expected3("1 2 3 1 1 3"), expected4("[(1, 2, 3); (1, 1, 3)]");
   threevector<int> v[2];
   v[0].set(1, 2, 3);
   v[1].set(1, 1, 3);

   outputVertex(v[0], str1);
   outputVertex(v[0], str2, "(", ")", ", ");
   ASSERT(expected1 == str1.str());
   ASSERT(expected2 == str2.str());

   outputVertices(&v[0], &v[2], str3);
   outputVertices(&v[0], &v[2], str4, "[", "]", "; ", "(", ")", ", ");
   ASSERT(expected3 == str3.str());
   ASSERT(expected4 == str4.str());

   for (int i = 0; i < clk.trefoilCoordsAsThreevectors.size(); i++)
      if(!testGenericOutput(clk.trefoilCoordsAsThreevectors[i])) return false;
      
   return true;
}

bool testRog(const vector<threevector<int> >& u, double expectedRog)
{
   double r2 = computeRog(u.begin(), u.end());
   double r = sqrt(r2);
   ASSERT_DOUBLES_EQUAL(expectedRog, r, DELTA);
   
   return true;
}

bool testRog()
{
	clktestclass clk;
 
   for (int i = 0; i < clk.trefoilCoordsAsThreevectors.size(); i++)
      if(!testRog(clk.trefoilCoordsAsThreevectors[i], clk.expectedRog[i])) return false;
   return true;
}

bool testWrAcn(const vector<threevector<int> >& u, double expectedWr, double expectedAcn)
{
   double Wr, Acn;
   computeWrAcn(u.begin(), u.end(), Wr, Acn);
   ASSERT_DOUBLES_EQUAL(expectedWr, Wr, DELTA);
   ASSERT_DOUBLES_EQUAL(expectedAcn, Acn, DELTA);
   
   return true;
}

bool testWrAcn()
{
	clktestclass clk;
 
   // First test that works correctly. No need to try on an empty range,
   // but it must work for one element ranges.
   ASSERT_EQUAL(0, findLast(0, 1));
   for (int i = -10; i < 10; i++)
      ASSERT_EQUAL(i, findLast(i, i + 1));
   // Now, a more thorough check of findLast.
   for (int i = -10; i < 10; i++)
      for (int j = i; j < i + 15; j++)
         ASSERT_EQUAL(j, findLast(i, j + 1));

   // "Finally ready to check Writhe and ACN.
   for (int i = 0; i < clk.trefoilCoordsAsThreevectors.size(); i++)
      if(!testWrAcn(clk.trefoilCoordsAsThreevectors[i], clk.expectedWr[i], clk.expectedAcn[i])) return false;
      
   return true;
}

bool testBfacf3()
{
   threevector<int> v;

   // Can we even instantiate clkConformationBfacf3 without crashing the system?
   clkCigar square;
   clkConformationBfacf3 knot(square);
   // Now, can we do anything with the conformation?
   ASSERT_EQUAL(1, knot.size());
   list<clkConformationAsList> components;
   knot.getComponents(components);
   ASSERT_EQUAL(1u, components.size());
   clkConformationAsList& comp0 = *components.begin();
   ASSERT_EQUAL(square.size(), comp0.size());
   ASSERT(comp0 == knot.getComponent(0));
   comp0.getVertex(0, v);
   v *= -1;
   comp0.translate(v);
   ASSERT(square == comp0);

   //   cout << "How about a two component link?" << endl;
   clkConformationAsList square2(square);
   square2.translate(0, 0, 1);
   clkConformationBfacf3 twoComponentLink(square, square2);
   ASSERT_EQUAL(2, twoComponentLink.size());
   list<clkConformationAsList> components2;
   twoComponentLink.getComponents(components2);
   ASSERT_EQUAL(2u, components2.size());
   list<clkConformationAsList>::iterator i = components2.begin();
   clkConformationAsList& comp1 = *i;
   i++;
   clkConformationAsList& comp2 = *i;
   ASSERT(comp1 == twoComponentLink.getComponent(0));
   ASSERT(comp2 == twoComponentLink.getComponent(1));
   comp1.getVertex(0, v);
   v *= -1;
   comp1.translate(v);
   comp2.translate(v);
   ASSERT(square == comp1);
   ASSERT(square2 == comp2);
   
   return true;
}

bool testBfacfProbabilitiesEqual(const bfacfProbabilities& p1, const bfacfProbabilities& p2)
{
   ASSERT_DOUBLES_EQUAL_MESSAGE("p0", p1.p0(), p2.p0(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_03p2", p1.p03p2(), p2.p03p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_2p0", p1.p2p0(), p2.p2p0(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_2p02p2", p1.p2p02p2(), p2.p2p02p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_4p2", p1.p4p2(), p2.p4p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_m23p2", p1.m23p2(), p2.m23p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_minus2", p1.m2(), p2.m2(), EPSILON);
   ASSERT_DOUBLES_EQUAL_MESSAGE("p_plus2", p1.p2(), p2.p2(), EPSILON);
   
   return true;
}

bool testBfacfProbabilitiesConsistent(const bfacfProbabilities& p)
{
   ASSERT_DOUBLES_EQUAL(4.0 * p.p2(), p.p4p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(p.p0() + 3.0 * p.p2(), p.p03p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(p.m2() + 3.0 * p.p2(), p.m23p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(2.0 * p.p0(), p.p2p0(), EPSILON);
   ASSERT_DOUBLES_EQUAL(2.0 * p.p0() + 2.0 * p.p2(), p.p2p02p2(), EPSILON);
   
   return true;
}

bool testBfacf3SetZ()
{
   // For first component.
   clkCigar square0;
   // For second component, translate square one unit in z-direction.
   clkConformationAsList square1(square0);
   square1.translate(0, 0, 1);

   double expectedDefaultZ = 0.20815;

   // Test Bfacf3 set z.

   // First try experiment a knot.
   clkConformationBfacf3 knot(square1);
   // The default z-value should be " << expectedDefaultZ << "." << endl;
   ASSERT_EQUAL(expectedDefaultZ, DEFAULT_Z);
   ASSERT_EQUAL(DEFAULT_Z, knot.getZ());
   ASSERT_EQUAL(DEFAULT_Z, knot.getComponent(0).getZ());

   // Check that the probabilities are what we expect.
   bfacfProbabilitiesFromZFixed pFromDefaultZ, p17(0.17), p19(0.19), p21(0.21);
   ASSERT_EQUAL(DEFAULT_Z, pFromDefaultZ.z());
   ASSERT_DOUBLES_EQUAL(0.88497198909891, pFromDefaultZ.m2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.03834267030036, pFromDefaultZ.p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.46165732969964, pFromDefaultZ.p0(), EPSILON);
   testBfacfProbabilitiesConsistent(pFromDefaultZ);

   ASSERT_DOUBLES_EQUAL(0.92021717125242, p17.m2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.02659427624919, p17.p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.47340572375081, p17.p0(), EPSILON);
   testBfacfProbabilitiesConsistent(p17);

   ASSERT_DOUBLES_EQUAL(0.90228277542182, p19.m2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.03257240819273, p19.p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.46742759180727, p19.p0(), EPSILON);
   testBfacfProbabilitiesConsistent(p19);

   ASSERT_DOUBLES_EQUAL(0.88315817362890, p21.m2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.03894727545703, p21.p2(), EPSILON);
   ASSERT_DOUBLES_EQUAL(0.46105272454297, p21.p0(), EPSILON);
   testBfacfProbabilitiesConsistent(p21);

   testBfacfProbabilitiesEqual(pFromDefaultZ, pFromDefaultZ);
   testBfacfProbabilitiesEqual(pFromDefaultZ, knot.getComponent(0).getProbabilities());

   // Try resetting z to 0.17 on the link level.
   knot.setZ(0.17);
   ASSERT_EQUAL(0.17, knot.getZ());
   ASSERT_EQUAL(0.17, knot.getComponent(0).getZ());
   testBfacfProbabilitiesEqual(p17, knot.getComponent(0).getProbabilities());

   // Try setting z to 0.21 on the component level.
   knot.getComponent(0).setZ(0.21);
   ASSERT_EQUAL(0.21, knot.getComponent(0).getZ());
   testBfacfProbabilitiesEqual(p21, knot.getComponent(0).getProbabilities());

   // Next try similar tests for a two-component link.
   clkConformationBfacf3 link(square0, square1);
   ASSERT_EQUAL(DEFAULT_Z, link.getZ());
   ASSERT_EQUAL(DEFAULT_Z, link.getComponent(0).getZ());
   ASSERT_EQUAL(DEFAULT_Z, link.getComponent(1).getZ());
   // Check probabilities for first component.
   testBfacfProbabilitiesEqual(pFromDefaultZ, link.getComponent(0).getProbabilities());
   // Check probabilities for second component.
   testBfacfProbabilitiesEqual(pFromDefaultZ, link.getComponent(1).getProbabilities());

   // Reset z to 0.17 on the link level.
   link.setZ(0.17);
   ASSERT_EQUAL(0.17, link.getZ());
   ASSERT_EQUAL(0.17, link.getComponent(0).getZ());
   ASSERT_EQUAL(0.17, link.getComponent(1).getZ());
   // Check probabilities for first component.
   testBfacfProbabilitiesEqual(p17, link.getComponent(0).getProbabilities());
   // Check probabilities for second component.
   testBfacfProbabilitiesEqual(p17, link.getComponent(1).getProbabilities());

   // Reset z to 0.19 on first component and 0.21 on second component.
   link.getComponent(0).setZ(0.19);
   link.getComponent(1).setZ(0.21);
   ASSERT_EQUAL(0.19, link.getComponent(0).getZ());
   ASSERT_EQUAL(0.21, link.getComponent(1).getZ());
   // Check probabilities for first component.
   testBfacfProbabilitiesEqual(p19, link.getComponent(0).getProbabilities());
   // Check probabilities for second component.
   testBfacfProbabilitiesEqual(p21, link.getComponent(1).getProbabilities());
   
   return true;
}

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

bool testRandomReset()
{
   int expected[] = {2, 46, 10, 80, 33, 41, 50, 1, 66, 80, 20, 71, 4, 41, 73, 71, 99, 70, 25, 4};
   list<int> l;
   // Can we save the state of our random number generator and restart later?


   // First verify that random number generator used by legacy is same as pseudorandom class.
   clkCigar square;
   clkConformationBfacf3 knot(square);
   knot.setSeed(42);
   pseudorandom r(42);
   for (int i = 0; i < 20; i++)
   {
      int m = rand_integer(1, 100);
      int n = r.rand_integer(1, 100);
      ASSERT_EQUAL(expected[i], m);
      ASSERT_EQUAL(m, n);
      l.push_back(m);
   }

   // Next verify that resetting both randomnumber generators with the same seed produces same result.
   r.sRandSimple(42);
   knot.setSeed(42);
   list<int>::const_iterator j = l.begin();
   for (int i = 0; i < 20; i++, j++)
   {
      int m = rand_integer(1, 100);
      int n = r.rand_integer(1, 100);
      ASSERT_EQUAL(m, n);
      ASSERT_EQUAL(m, *j);
   }

   // Now, try saving both states, run different sequences and see if we can restart.

   // save states
   pseudorandom s(r);
   pseudorandom t;
   t = saveRandomState();

   // save what would be results running from this state.
   l.clear();
   for (int i = 0; i < 20; i++)
      l.push_back(r.rand_integer(1, 100));

   // scramble random number generators.
   r.sRandSimple(907);
   knot.setSeed(481);
   for (int i = 0; i < 1000; i++)
   {
      r.rand_integer(1, 100);
      rand_integer(1, 100);
   }

   copyRandomState(t);
   j = l.begin();
   for (int i = 0; i < 20; i++, j++)
   {
      int m = s.rand_integer(1, 100);
      int n = t.rand_integer(1, 100);
      int o = rand_integer(1, 100);
      ASSERT_EQUAL(m, n);
      ASSERT_EQUAL(m, o);
      ASSERT_EQUAL(m, *j);
   }
   
   return true;
}

bool testBfacf3Run()
{
   threevector<int> v;

   // The crucial test is whether we can actually run the BFACF algorithm in a repeatable way.
   string targetUnkot("0 1 0 -1 1 0 -1 1 -1 -1 2 -1 -1 3 -1 0 3 -1 0 3 0 0 2 0 -1 2 0 -2 2 0 -2 1 0 -3 1 0 -3 0 0 -2 0 0 -1 0 0 0 0 0");
   clkCigar square;
   clkConformationAsList targetConformation;
   targetConformation.readFromText(targetUnkot);
   targetConformation.getVertex(0, v);
   v *= -1;
   targetConformation.translate(v);
   // Set seed to 42 and apply 100 BFACF moves. The resulting conformation should be targetConformation.

   // Instantiate two bfacf classes with same initial conformations.
   clkConformationBfacf3 knot0(square), knot1(square);
   // Set seed for first conformation to 42 and perform 100 steps.
   knot0.setSeed(42);
   for (int i = 0; i < 100; i++)
      knot0.step();
   clkConformationAsList saved(knot0.getComponent(0));
   clkConformationAsList result0(knot0.getComponent(0));
   result0.getVertex(0, v);
   v *= -1;
   result0.translate(v);

   // Set seed for second conformation to 42 and perform 100 steps.
   knot1.setSeed(42);
   for (int i = 0; i < 100; i++)
      knot1.step();
   clkConformationAsList result1(knot1.getComponent(0));
   result1.getVertex(0, v);
   v *= -1;
   result1.translate(v);

   ASSERT_MESSAGE("knot0 and knot1 should be same.", knot0.getComponent(0) == knot1.getComponent(0));
   ASSERT_MESSAGE("knot0 should not have when iterating knot1.", knot0.getComponent(0) == saved);
   ASSERT_MESSAGE("Both knots should be same as target conformation.", result0 == targetConformation);
   ASSERT(result0 == result1);

   // Can we save the state of the bfacf3 and restart at this same point?
   // run another 1000 steps
   clkConformationBfacf3 knot2(targetConformation);
   pseudorandom saveRandom = saveRandomState();
   for (int i = 0; i < 100; i++)
      knot2.step();
   clkConformationAsList secondTargetCoformation(knot2.getComponent(0));
   // try restoring random number generator and run bfacf again
   clkConformationBfacf3 knot3(targetConformation);
   copyRandomState(saveRandom);
   for (int i = 0; i < 100; i++)
      knot3.step();
   clkConformationAsList actualConformation(knot3.getComponent(0));
   secondTargetCoformation.getVertex(0, v);
   v *= -1;
   secondTargetCoformation.translate(v);
   actualConformation.getVertex(0, v);
   v *= -1;
   actualConformation.translate(v);
   ASSERT(secondTargetCoformation == actualConformation);
   
   return true;
}

bool testReadFromText()
{
	clktestclass clk;
 
   stringstream stream;
   stream << clk.manyTrefoilString;
   clkConformationAsList knot;

   // we should be able to read exactly five conformations from this stream
   ASSERT_MESSAGE("Could not read trefoil0", knot.readFromText(stream));
   ASSERT_MESSAGE("trefoil0 does not match expected conformation.", knot == *clk.clkTrefoil0);
   
   ASSERT_MESSAGE("Could not read trefoil1", knot.readFromText(stream));
   ASSERT_MESSAGE("trefoil1 does not match expected conformation.", knot == *clk.clkTrefoil1);
   
   ASSERT_MESSAGE("Could not read trefoil2", knot.readFromText(stream));
   ASSERT_MESSAGE("trefoil2 does not match expected conformation.", knot == *clk.clkTrefoil2);
   
   ASSERT_MESSAGE("Could not read trefoil3", knot.readFromText(stream));
   ASSERT_MESSAGE("trefoil3 does not match expected conformation.", knot == *clk.clkTrefoil3);
   
   ASSERT_MESSAGE("Could not read trefoil4", knot.readFromText(stream));
   ASSERT_MESSAGE("trefoil4 does not match expected conformation.", knot == *clk.clkTrefoil4);
   
   ASSERT_MESSAGE("There should only be five conformations to read.", !knot.readFromText(stream));
   
   return true;
}

bool testReadFromCoords()
{
	clktestclass clk;

   stringstream stream;
   stream << clk.manyTrefoilCoordString;
   clkConformationAsList knot;

   // we should be able to read exactly five conformations from this stream
   ASSERT_MESSAGE("Could not read trefoil0", knot.readFromCoords(stream));
   ASSERT_MESSAGE("trefoil0 does not match expected conformation.", knot == *clk.clkTrefoil0);
   
   ASSERT_MESSAGE("Could not read trefoil1", knot.readFromCoords(stream));
   ASSERT_MESSAGE("trefoil1 does not match expected conformation.", knot == *clk.clkTrefoil1);
   
   ASSERT_MESSAGE("Could not read trefoil2", knot.readFromCoords(stream));
   ASSERT_MESSAGE("trefoil2 does not match expected conformation.", knot == *clk.clkTrefoil2);
   
   ASSERT_MESSAGE("Could not read trefoil3", knot.readFromCoords(stream));
   ASSERT_MESSAGE("trefoil3 does not match expected conformation.", knot == *clk.clkTrefoil3);
   
   ASSERT_MESSAGE("Could not read trefoil4", knot.readFromCoords(stream));
   ASSERT_MESSAGE("trefoil4 does not match expected conformation.", knot == *clk.clkTrefoil4);
   
   ASSERT_MESSAGE("There should only be five conformations to read.", !knot.readFromText(stream));
   
   return true;
}

bool testBfacf3CountEdges()
{
   clkCigar square;
   clkConformationAsList square2(square);
   square2.translate(0,0,1);
   clkConformationBfacf3 onecomp(square);
   //ASSERT_EQUAL(4, onecomp.totalEdges());
   ASSERT(false);
   clkConformationBfacf3 twocomps(square, square2);
   //ASSERT_EQUAL(8, twocomps.totalEdges());
   ASSERT(false);
   
   return true;
}

bool testBfacf3WithQ(){
	//need to test BFACF3 with Q =/= 1 here.
	clkCigar square;
	clkConformationAsList square1(square);
	clkConformationBfacf3 knot0(square1);
	clkConformationBfacf3 knot1(square1);

	ASSERT(knot0.size() == knot1.size());
	knot0.setZ(.2000);
	knot1.setZ(.2000);
	ASSERT(knot0.getZ() == knot1.getZ());
	knot0.step(12345);
	knot1.stepQ(12345, 1, .2100);
	//would prefer to compare the actual vectors of each conformation. 
	ASSERT(knot0.size() == knot1.size());
	
	return true;
}

bool testPrecomputedBfacf3Probs(){
	clkCigar square;
	clkConformationAsList square1(square);
	testClkConformationBfacf3 knot(square1);
	//select a z with q=1
	double z = .2000, q = 1;
	//precompute probabilities
	knot.init_Q(z, q);
	//compare precomputed probabilities with directly computed probabilities
	for (int n = 4; n < MAX_PRECOMPUTE_LENGTH; n++){
		ASSERT(knot.probMap[n].p_plus2 == (pow((n + 2), (q - 1))*(z * z)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z));
		ASSERT(knot.probMap[n].p_minus2 == pow(n, (q - 1)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z));
		ASSERT(knot.probMap[n].p_0 == .5*((pow((n + 2), (q - 1))*(z * z)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z) +
			pow(n, (q - 1)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z)));
	}
	//select a different z with q=3
	z = .1812; q = 3;
	//precompute
	knot.init_Q(z, q);
	//compare
	for (int n = 4; n < MAX_PRECOMPUTE_LENGTH; n++){
		ASSERT(knot.probMap[n].p_plus2 == (pow((n + 2), (q - 1))*(z * z)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z));
		ASSERT(knot.probMap[n].p_minus2 == pow(n, (q - 1)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z));
		ASSERT(knot.probMap[n].p_0 == .5*((pow((n + 2), (q - 1))*(z * z)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z) +
			pow(n, (q - 1)) / (pow(n, (q - 1)) + 3.0*pow((n + 2), q - 1) * z * z)));
	}
	
	return true;
}

