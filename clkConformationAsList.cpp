/* 
 * File:   clkConformationAsList.cpp
 * Author: kmo
 * 
 * Created on December 4, 2012, 3:04 PM
 */

#include "clkConformationAsList.h"

#include "genericConformation.h"

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

clkConformationAsList::clkConformationAsList() { }

clkConformationAsList::clkConformationAsList(const clk& orig)
{
   threevector<int> v;
   int n = orig.size();
   for (int i = 0; i < n; i++)
   {
      orig.getVertex(i, v);
      data.push_back(v);
   }
}

clkConformationAsList::clkConformationAsList(const int* coords, int numVerts)
{
   threevector<int> v;
   for (int i = 0; i < numVerts; i++)
   {
      v.set(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
      data.push_back(v);
   }
}

clkConformationAsList::clkConformationAsList(const clkConformationAsList& orig) : data(orig.data) { }

clkConformationAsList::~clkConformationAsList() { }

int clkConformationAsList::size() const
{
   return data.size();
}

// this function is not efficient

void clkConformationAsList::getVertex(int index, threevector<int>& v) const
{
   deque<threevector<int> >::const_iterator iterator = data.begin();
   for (int i = 0; i < index; i++)
      iterator++;
   v = *iterator;
}

void clkConformationAsList::addVertexBack(const threevector<int>& v)
{
   data.push_back(v);
}

void clkConformationAsList::addVertexBack(int x, int y, int z)
{
   threevector<int> v(x, y, z);
   addVertexBack(v);
}

void clkConformationAsList::addVertexFront(const threevector<int>& v)
{
   data.push_front(v);
}

void clkConformationAsList::addVertexFront(int x, int y, int z)
{
   threevector<int> v(x, y, z);
   addVertexFront(v);
}

void clkConformationAsList::clear()
{
   data.clear();
}

void clkConformationAsList::writeAsText(ostream& os) const
{
   deque<threevector<int> >::const_iterator i = data.begin();
   while (i != data.end())
   {
      if (i != data.begin()) os << " ";
      os << *i;
      i++;
   }
   os << endl;
}

string clkConformationAsList::writeAsText() const
{
   stringstream str;
   writeAsText(str);
   return str.str();
}

bool clkConformationAsList::readFromText(istream& is)
{
   clear();

   string line;
   deque<int> numbers;

   while (numbers.size() == 0)
   {
      if (!is) return false;
      getline(is, line);
      stringstream str;
      str << line;
      int n;
      while (str)
      {
         if (str >> n)
            numbers.push_back(n);
      }
   }

   int n = numbers.size();
   if (n < 3 || n % 3 == 2)
      return false;

   deque<int>::iterator i = numbers.begin();
   if (numbers.size() % 3 == 1)
   {
      // numbers is a list of integers whose size is a positive multiple of three
      // plus one or two
      // assume that the first entry is the number of vertices and ignore it
      i++;
   }

   while (i != numbers.end())
   {
      int x = *i;
      i++;
      int y = *i;
      i++;
      int z = *i;
      i++;
      addVertexBack(x, y, z);
   }

   return true;
}

bool clkConformationAsList::readFromText(const string& text)
{
   stringstream str;
   str << text;
   return readFromText(str);
}

bool clkConformationAsList::lineIsVertex(std::string line)
{
   int count = 0;
   stringstream str;
   str << line;
   int n;
   while (str)
   {
      if (str >> n)
         count++;
   }

   return count == 3;
}

bool clkConformationAsList::readFromCoords(std::istream& is)
{
   clear();
   if (!is) return false;

   string line;
   stringstream singleLine;

   getline(is, line);
   while (lineIsVertex(line))
   {
      singleLine << line << " ";
      getline(is, line);
      if (!is) break;
   }

   //   cout << singleLine.str() << endl;

   return readFromText(singleLine);
}


// to mimic BitFile structure and related functions

#define MOVE_NORTH 0
#define MOVE_EAST  1
#define MOVE_WEST  2
#define MOVE_SOUTH 3
#define MOVE_UP    4
#define MOVE_DOWN  5

#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

class bitio
{
public:

   bitio() : ch(0), index(0) { }

   int read_bit(istream& is)
   {
      if (index % 8 == 0)
      {
         is.read((char*) &ch, 1);
         index = 0;
      }

      int value = ch >> index;
      index++;

      return 0x1 & value;

   }

   void write_bit(int bit, ostream& os)
   {
      ch |= (bit & 0x1) << index;
      if (++index == 8)
      {
         os.write((char*) &ch, 1);
         ch = 0;
         index = 0;
      }

   }

   void write_increment(int value, ostream& os)
   {
      write_bit(value, os);
      write_bit(value >> 1, os);
      write_bit(value >> 2, os);
   }

   int read_increment(istream& is)
   {
      int value = 0;
      for (int i = 0; i < 3; i++)
         value |= read_bit(is) << i;

      return value;
   }

   bool read_increment(threevector<int>& incr, istream& is)
   {
      int dir;
      switch (dir = read_increment(is))
      {
         case MOVE_NORTH:
            incr.set(0, 1, 0);
            break;
         case MOVE_EAST:
            incr.set(1, 0, 0);
            break;
         case MOVE_WEST:
            incr.set(-1, 0, 0);
            break;
         case MOVE_SOUTH:
            incr.set(0, -1, 0);
            break;
         case MOVE_DOWN:
            incr.set(0, 0, -1);
            break;
         case MOVE_UP:
            incr.set(0, 0, 1);
            break;
         default:
            cerr << "bad move direction (" << dir << ") found!" << endl;
            return false;
      }
      //  fprintf (stderr, "%c", dirname [dir]);

      return true;
   }

   bool write_increment(const threevector<int>& incr, ostream& os)
   {
      int value;
      switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
      {
         case NORTH:
            value = MOVE_NORTH;
            break;
         case EAST:
            value = MOVE_EAST;
            break;
         case WEST:
            value = MOVE_WEST;
            break;
         case SOUTH:
            value = MOVE_SOUTH;
            break;
         case UP:
            value = MOVE_UP;
            break;
         case DOWN:
            value = MOVE_DOWN;
            break;
         default:
            cerr << "bad increment (" <<
                    ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
                    << ") (" << incr << ")" << endl;
            return false;
      }
      write_increment(value, os);
      //  fprintf (stderr, "%c", dirname [value]);
      return true;
   }

private:
   unsigned char ch;
   int index;

};

class intio
{
public:

   intio() { }

   void write(int val, ostream& os)
   {
      unsigned char a, b, c, d;

      a = (unsigned char) (0xff & val >> 24);
      b = (unsigned char) (0xff & val >> 16);
      c = (unsigned char) (0xff & val >> 8);
      d = (unsigned char) (0xff & val);

      os.put(a);
      os.put(b);
      os.put(c);
      os.put(d);
   }

   void write(const threevector<int>& v, ostream& os)
   {
      write(v.getX(), os);
      write(v.getY(), os);
      write(v.getZ(), os);
   }

   bool read(int& v, istream& is)
   {
      if (!is) return false;
      unsigned char uc;
      v = 0;
      for (int shift = 24; shift >= 0; shift -= 8)
      {
         is.get((char&) uc);
         v |= uc << shift;
      }
      return true;
   }

   bool read(threevector<int>& v, istream& is)
   {
      int x, y, z;
      if (!read(x, is)) return false;
      if (!read(y, is)) return false;
      if (!read(z, is)) return false;

      v.set(x, y, z);
      return true;
   }

};

bool clkConformationAsList::writeAsCube(std::ostream& os) const
{
   // code in this function modified by Reuben Brasher 
   // from Rob Sharein's seqconvert 
   intio intio;

   int num_vertices = size();
   int numbits = num_vertices * 3;
   int numbytes = (numbits + 7) / 8;

   os.write("CUBE", 4);
   int S; // size of field is 16 bytes (number of vertices plus starting point) plus numbytes
   S = 16 + numbytes;
   intio.write(S, os);
   //   os.write((const char*) &S, sizeof (int));

   // write number of vertices
   intio.write(num_vertices, os);
   //   os.write((const char*) &num_vertices, sizeof (int));

   // write starting location
   deque<threevector<int> >::const_iterator i = data.begin();
   //   int firstVertex[3] = {i->getX(), i->getY(), i->getZ()};
   //   os.write((const char*) &firstVertex[0], 3 * sizeof (int));
   intio.write(*i, os);

   bitio bf;

   while (i != data.end())
   {
      threevector<int> increment(*i);
      i++;
      if (i != data.end())
         increment -= *i;
      else
         increment -= *data.begin();
      increment *= -1;
      if (!bf.write_increment(increment, os)) return false;
   }

   int num_pad_bits = numbytes * 8 - numbits;

   for (int pad = 0; pad < num_pad_bits; pad++)
      bf.write_bit(1, os);

   return true;

}

bool clkConformationAsList::readFromCube(std::istream& is)
{
   // code in this function modified by Reuben Brasher 
   // from Rob Sharein's seqconvert 
   char tag_name [5];
   intio intio;
   if (!is.read(tag_name, 4)) return false;
   tag_name [4] = (char) 0;
   string cube("CUBE");

   if (cube != tag_name)
      cerr << "bad tag >" << tag_name << "<" << endl;

   int S, num_vertices;
   threevector<int> v;

   intio.read(S, is);
   intio.read(num_vertices, is);
   intio.read(v, is);

   clear();
   addVertexBack(v);

   bitio bf;
   threevector<int> increment;

   for (int i = 1; i < num_vertices; i++)
   {
      if (!bf.read_increment(increment, is)) return false;
      v += increment;
      addVertexBack(v);
   }
   if (!bf.read_increment(increment, is)) return false; // ignore 

   int numbits = num_vertices * 3;
   int numbytes = (numbits + 7) / 8;
   int num_pad_bits = numbytes * 8 - numbits;

   for (int pad = 0; pad < num_pad_bits; pad++)
   {
      if (bf.read_bit(is) != 1)
      {
         cerr << "bad bit seen!" << endl;
         //      fflush (stderr);
         return false;
      }
   }

   return true;
}

void clkConformationAsList::translate(const threevector<int>& offset)
{
   deque<threevector<int> >::iterator i = data.begin();
   while (i != data.end())
   {
      *i += offset;
      i++;
   }
}

void clkConformationAsList::translate(int x, int y, int z)
{
   translate(threevector<int>(x, y, z));
}

double clkConformationAsList::rog() const
{
   return computeRog(data.begin(), data.end());
}

void clkConformationAsList::writheACN(double& writhe, double& acn) const
{
   computeWrAcn(data.begin(), data.end(), writhe, acn);
}

void clkConformationAsList::writheACN(const clkConformationAsList& other, double& writhe, double& acn) const
{
   computeWrAcnTwocomps(data.begin(), data.end(), other.data.begin(), other.data.end(), writhe, acn);
}