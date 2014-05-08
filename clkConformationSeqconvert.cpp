/* 
 * File:   clkConformationSeqconvert.cpp
 * Author: kmo
 * 
 * Created on December 4, 2012, 3:08 PM
 */

#include <string.h>

#include "clkConformationSeqconvert.h"
#include "legacy.h"

using namespace std;

class clkConformationSeqconvert::impl
{

   impl() : coord(NULL), length(0), maxlength(0) { }

   impl(const clk& orig) : length(orig.size()), maxlength(orig.size())
   {
      coord = new ivector[length];
      //   data = coord;
      threevector<int> u;
      for (int i = 0; i < length; i++)
      {
         orig.getVertex(i, u);
         set_ivector(coord[i], u.getX(), u.getY(), u.getZ());
      }
   }

   impl(const clkConformationSeqconvert& orig) : length(orig.size()), maxlength(orig.size())
   {
      ivector* coord = new ivector[length];
      memcpy(coord, orig.implementation->coord, sizeof (ivector) * length);
   }

   ~impl()
   {
      if (coord) delete[] coord;
   }

   void copy(const clk& orig)
   {
      if (coord) delete[] coord;

      length = maxlength = orig.size();
      //   maxlength = length;
      coord = new ivector[length];
      threevector<int> u;
      for (int i = 0; i < length; i++)
      {
         orig.getVertex(i, u);
         set_ivector(coord[i], u.getX(), u.getY(), u.getZ());
      }
   }

   void copy(const clkConformationSeqconvert& orig)
   {
      if (coord) delete[] coord;

      length = maxlength = orig.size();
      //   maxlength = length;
      coord = new ivector[length];
      memcpy(coord, orig.implementation->coord, sizeof (ivector) * length);
   }

   ivector* coord;
   int length;
   int maxlength;

   friend class clkConformationSeqconvert;
};

clkConformationSeqconvert::clkConformationSeqconvert():
implementation(new clkConformationSeqconvert::impl) { }

clkConformationSeqconvert::clkConformationSeqconvert(const clk& orig) :
implementation(new clkConformationSeqconvert::impl(orig)) { }

clkConformationSeqconvert::clkConformationSeqconvert(const clkConformationSeqconvert& orig) :
implementation(new clkConformationSeqconvert::impl(orig)) { }

clkConformationSeqconvert::~clkConformationSeqconvert()
{
   delete implementation;
}

int clkConformationSeqconvert::size() const
{
   return implementation->length;
}

void clkConformationSeqconvert::getVertex(int index, threevector<int>& v) const
{
   v.set(implementation->coord[index][0],
           implementation->coord[index][1],
           implementation->coord[index][2]);
}

clkConformationSeqconvert& clkConformationSeqconvert::operator=(const clk& orig)
{
   implementation->copy(orig);

   return *this;
}

clkConformationSeqconvert& clkConformationSeqconvert::operator=(const clkConformationSeqconvert& orig)
{
   implementation->copy(orig);

   return *this;
}

void clkConformationSeqconvert::computeWritheACN(double& writhe, double& ACN)
{
   writhe = ::writhe(ACN, size(), implementation->coord, 0);
}

double clkConformationSeqconvert::radiusOfGyration()
{
   return radius_of_gyration(size(), implementation->coord);
}

