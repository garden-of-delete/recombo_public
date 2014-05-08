/* 
 * File:   mainfunctions.cpp
 * Author: kmo
 *
 * Created on March 16, 2013, 3:37 PM
 */

#include <iostream>
#include <fstream>

#include "recomboArgprocessor.h"
#include "runRecombo.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "random.h"
#include "autocorr.h"
#include "conformationAsList.h"

using namespace std;

int mainRecombo(int argc, char** argv)
{
   //   cout << "Hello, recombo." << endl;
   recomboArgprocessor argprocessor;
   if (!argprocessor.process(argc, argv))
   {
      cout << argprocessor;
      return 1;
   }

   runRecombo *run;
   if (argprocessor.aProtocol())
      run = new runRecomboA(argprocessor);
   else if (argprocessor.hProtocol())
      run = new runRecomboH(argprocessor);
   else
      run = new runRecomboB(argprocessor);

   if (!run->readInitialCoords())
   {
      delete run;
      return 2;
   }

   run->writeLog(cout);

   if (!run->openfiles())
   {
      //      cerr << "Could not open files for output." << endl;
      delete run;
      return 3;
   }

   run->dorecombo();

   run->closefiles();

   delete run;
   return 0;
}

int mainRerecombo(int argc, char** argv)
{
   cout << "Hello, recombo." << endl;
   rerecomboArgprocessor argprocessor;
   if (!argprocessor.process(argc, argv))
   {
      cout << argprocessor;
      return 1;
   }

   runRerecombo run(argprocessor);

   run.writeLog(cout);

   if (!run.openfiles())
   {
//      cerr << "Could not open files for output." << endl;
      return 3;
   }

   run.dorecombo();

   run.closefiles();

   return 0;
}

int mainTwoCompWrithe(int argc, char** argv)
{
   if (argc < 2)
   {
      cout << "Usage: twocompwrithe inputfilename" << endl;
      return 1;
   }

   fstream data(argv[1], ios::in | ios::binary);
   if (!data)
   {
      cerr << "Could not open " << argv[1] << " for input." << endl;
      return 2;
   }

//   ostream *stats;
//   if (argc > 2)
//   {
//      stats = new ofstream(argv[2]);
//      if (!(*stats))
//      {
//         cerr << "Could not open " << argv[2] << " for output." << endl;
//         data.close();
//         delete stats;
//         return 3;
//      }
//      else
//      {
//         stats = &cout;
//      }
//   }

   double wr0, wr1, wr01, acn;
   vector<double> d0, d1, d01, dtotal, dcross;
   autocorr ac;
   double m0, v0, m1, v1, m01, v01, mtotal, vtotal, covar;

   clkConformationAsList comp0, comp1;

   int count = -1;

   //   cout << "Columns contain writhe for component 0, component 1 and the writhe of the link in that order." << endl;

   while (count && data)
   {
//      if (!comp0.readFromText(data))
      if (!comp0.readFromCube(data))
      {
         //         cerr << "Error reading comp0." << endl;
         break;
      }
      if (!data)
      {
         cerr << "Unexpected end of file." << endl;
         break;
      }
//      if (!comp1.readFromText(data))
      if (!comp1.readFromCube(data))
      {
         cerr << "Error reading comp1." << endl;
         break;
      }
//      *stats << "comp0 " << comp0 << endl << "comp1 " << comp1 << endl;
      comp0.writheACN(wr0, acn);
//      *stats << wr0 << "\t";
      comp1.writheACN(wr1, acn);
//      *stats << wr1 << "\t";
      comp0.writheACN(comp1, wr01, acn);
//      *stats << wr01 << "\t" << wr0 + wr1 + wr01;
//      *stats << endl;
      d0.push_back(wr0);
      d1.push_back(wr1);
      d01.push_back(wr01);
      dtotal.push_back(wr0 + wr1 + wr01);
      dcross.push_back(wr0 * wr1);

      count--;
   }
   data.close();
//   if(stats && stats != &cout)
//   {
//      stats->close();
//      delete stats;
//   }

   m0 = ac.computeMean(d0);
   v0 = ac.computeVariance(d0);
   m1 = ac.computeMean(d1);
   v1 = ac.computeVariance(d1);
   m01 = ac.computeMean(d01);
   v01 = ac.computeVariance(d01);
   mtotal = ac.computeMean(dtotal);
   vtotal = ac.computeVariance(dtotal);
   covar = ac.computeMean(dcross) - m0 * m1;
   //   cout << "Mean and standard deviation of writhe for component 0:" << endl;
   //   cout << m0 << "\t" << std0 << endl;
   //   cout << "Mean and standard deviation of writhe for component 1:" << endl;
   //   cout << m1 << "\t" << std1 << endl;
   //   cout << "Mean and standard deviation of writhe between components:" << endl;
   //   cout << m01 << "\t" << std01 << endl;
   //   cout << "Mean and standard deviation of writhe for total writhe:" << endl;
   //   cout << mtotal << "\t" << stdtotal << endl;
   cout << m0 << "\t" << v0 << "\t"
           << m1 << "\t" << v1 << "\t"
           << covar << "\t"
           << m01 << "\t" << v01 << "\t"
           << mtotal << "\t" << vtotal << endl;
   return 0;
}

int mainTestCube(int arc, char** argv)
{
   if (arc < 2)
      return 1;
   fstream cubes(argv[1]);
   if (!cubes)
   {
      cerr << "Could not open " << argv[1] << " for input." << endl;
      return 1;
   }

   clkConformationAsList comp;
   int count = 0;
   while (cubes)
   {
      if (!comp.readFromCube(cubes))
      {
         cerr << "asdf" << endl;
         break;
      }
      cout << ++count << endl << comp << endl;
   }
   cubes.close();
   return 0;
}

int mainFixegc(int argc, char**argv)
{
   if (argc < 2)
   {
      cerr << "Usage: fixegc filename.egc" << endl;
      return 1;
   }

   ifstream in(argv[1]);
   while (in)
   {
      string line;
      if (getline(in, line))
      {
         if (line == "")
            cout << "b1-a1-" << endl;
         else if (line == "|")
            cout << "a1+a2-|b1+b2-" << endl;
         else cout << line << endl;
      }
   }
   in.close();

   return 0;
}

int mainWrithe(int argc, char**argv)
{
   if (argc < 2)
   {
      cerr << "Usage: writhe filename.txt" << endl;
      return 1;
   }

   ifstream in(argv[1]);
   if(!in)
   {
      cerr << "Could not open " << argv[1] << " for input." << endl;
      return 2;
   }
   
   conformationAsList con;
   double wr, acn;
   while (in)
   {
      if(con.readFromText(in))
      {
//         cout << con << endl;
         con.writheACN(wr, acn);
         cout << wr << endl;
      }
   }
   
   return 0;
}