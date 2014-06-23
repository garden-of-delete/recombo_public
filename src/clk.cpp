/* 
 * File:   clk.cpp
 * Author: kmo
 * 
 * Created on November 27, 2012, 6:49 PM
 */

#include "clk.h"

using namespace std;

clk::clk()
{ }

clk::clk(const clk& orig): conformation(orig)
{ }

clk::~clk() { }

void clk::getVertex(int index, threevector<double>& v) const
{
   threevector<int> vi;
   getVertex(index, vi);
   v.set(vi.getX(), vi.getY(), vi.getZ());
}

void clk::display(ostream& out) const
{
   int n = size();
   threevector<int> v;

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

void clk::getVertices(list<threevector<int> >& list) const
{
   int n = size();
   for (int i = 0; i < n; i++)
   {
      threevector<int> v;
      getVertex(i, v);
      list.push_back(v);
   }

}

bool clk::operator==(const clk& other) const
{
   int n = size();
   if(n != other.size()) return false;
   threevector<int> u, v;
   for(int i = 0; i < n; i++)
   {
      getVertex(i, u);
      other.getVertex(i, v);
      if(u != v) return false;
   }
   return true;
}

bool clk::operator !=(const clk& other) const
{
   return !(*this == other);
}