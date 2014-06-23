/* 
 * File:   conformation.cpp
 * Author: kmo
 * 
 * Created on November 30, 2012, 6:49 PM
 */

#include "conformation.h"

using namespace std;

conformation::conformation() { }

conformation::conformation(const conformation& orig) { }

conformation::~conformation() { }

void conformation::display(ostream& out) const
{
   int n = size();
   threevector<double> v;

   if (n > 0)
   {
      getVertex(0, v);
      out << v;
   }
   for (int i = 1; i < n; i++)
   {
      getVertex(i, v);
      out << ' ' << v;
   }
}

void conformation::getVertices(list<threevector<double> >& list) const
{
   int n = size();
   for (int i = 0; i < size(); i++)
   {
      threevector<double> v;
      getVertex(i, v);
      list.push_back(v);
   }
}

ostream& operator<<(ostream& out, const conformation& conformation)
{
   conformation.display(out);
   return out;
}
