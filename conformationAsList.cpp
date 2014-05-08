/* 
 * File:   conformationAsList.cpp
 * Author: kmo
 * 
 * Created on December 4, 2012, 3:03 PM
 */

#include "conformationAsList.h"
#include "genericConformation.h"

#include <deque>
#include <sstream>

using namespace std;

conformationAsList::conformationAsList() { }

conformationAsList::conformationAsList(const conformationAsList& orig) : data(orig.data) { }

conformationAsList::conformationAsList(const conformation& orig)
{
   threevector<double> v;
   int n = orig.size();
   for (int i = 0; i < n; i++)
   {
      orig.getVertex(i, v);
      data.push_back(v);
   }
}

conformationAsList::~conformationAsList() { }

int conformationAsList::size() const
{
   return data.size();
}

// this function is not efficient

void conformationAsList::getVertex(int index, threevector<double>& v) const
{
   vector<threevector<double> >::const_iterator iterator = data.begin();
   for (int i = 0; i < index; i++)
      iterator++;
   v = *iterator;
}

void conformationAsList::addVertexBack(const threevector<double>& v)
{
   data.push_back(v);
}

void conformationAsList::addVertexBack(double x, double y, double z)
{
   threevector<double> v(x, y, z);
   addVertexBack(v);
}

void conformationAsList::clear()
{
   data.clear();
}

void conformationAsList::writeAsText(ostream& os) const
{
   vector<threevector<double> >::const_iterator i = data.begin();
   while (i != data.end())
   {
      if (i != data.begin()) os << " ";
      os << *i;
      i++;
   }
   os << endl;
}

string conformationAsList::writeAsText() const
{
   stringstream str;
   writeAsText(str);
   return str.str();
}

bool conformationAsList::readFromText(istream& is)
{
   clear();

   string line;
   deque<double> numbers;

   while (numbers.size() == 0)
   {
      if (!is) return false;
      getline(is, line);
      stringstream str;
      str << line;
      double n;
      while (str)
      {
         if (str >> n)
            numbers.push_back(n);
      }
   }

   int n = numbers.size();
   if (n < 3 || n % 3 == 2)
      return false;

   deque<double>::iterator i = numbers.begin();
   if (numbers.size() % 3 == 1)
   {
      // numbers is a list of integers whose size is a positive multiple of three
      // plus one or two
      // assume that the first entry is the number of vertices and ignore it
      i++;
   }

   while (i != numbers.end())
   {
      double x = *i;
      i++;
      double y = *i;
      i++;
      double z = *i;
      i++;
      addVertexBack(x, y, z);
   }

   return true;
}

bool conformationAsList::readFromText(const string& text)
{
   stringstream str;
   str << text;
   return readFromText(str);
}

void conformationAsList::translate(const threevector<double>& offset)
{
   vector<threevector<double> >::iterator i = data.begin();
   while (i != data.end())
   {
      *i += offset;
      i++;
   }
}

void conformationAsList::translate(double x, double y, double z)
{
   translate(threevector<double>(x, y, z));
}

double conformationAsList::rog() const
{
   return computeRog(data.begin(), data.end());
}

void conformationAsList::writheACN(double& writhe, double& acn) const
{
   computeWrAcn(data.begin(), data.end(), writhe, acn);
}
