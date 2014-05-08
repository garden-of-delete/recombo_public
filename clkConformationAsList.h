/* 
 * File:   clkConformationAsList.h
 * Author: kmo
 *
 * Created on December 4, 2012, 3:04 PM
 */

#ifndef CLKCONFORMATIONASLIST_H
#define	CLKCONFORMATIONASLIST_H

#include "clk.h"

//#include <list>
#include <deque>
#include <string>

/**
 * Conformation implemented using a container of vertices. Vertices can be added
 * removed or changed.
 */
class clkConformationAsList : public clk
{
public:
   /**
    * Default constructor. After construction, vertices should be added before
    * conformation accessed. Otherwise, there is no guarantee that conformation
    * will have any particular geometry.
    */
   clkConformationAsList();

   /**
    * Constructor that reads vertices from an array of integers.
    * @param coords an array of integers assumed to be of length 3 * numVerts
    * @param numVerts number of vertices
    */
   clkConformationAsList(const int* coords, int numVerts);

   /**
    * Copy constructor.
    * @param orig new conformation will have same coordinates as orig.
    */
   clkConformationAsList(const clk& orig);

   /**
    * Copy constructor.
    * @param orig new conformation will have same coordinates as orig.
    */
   clkConformationAsList(const clkConformationAsList& orig);

   /**
    * Destructor.
    */
   virtual ~clkConformationAsList();

   virtual int size() const;
   virtual void getVertex(int index, threevector<int>& v) const;

   /**
    * Add a new vertex in new last position of conformation. This operation will
    * increase size by one.
    * @param v vertex to add.
    */
   virtual void addVertexBack(const threevector<int>& v);

   /**
    * Add a new vertex in new last position of conformation. This operation will
    * increase size by one.
    * @param x first coordinate of vertex to add.
    * @param y second coordinate of vertex to add.
    * @param z third coordinate of vertex to add.
    */
   virtual void addVertexBack(int x, int y, int z);

   /**
    * Add a new vertex in new first position of conformation. This operation 
    * will increase size by one and shift position of all existing vertices by
    * one.
    * @param v vertex to add.
    */
   virtual void addVertexFront(const threevector<int>& v);

   /**
    * Add a new vertex in new first position of conformation. This operation 
    * will increase size by one and shift position of all existing vertices by
    * one.
    * @param x first coordinate of vertex to add.
    * @param y second coordinate of vertex to add.
    * @param z third coordinate of vertex to add.
    */
   virtual void addVertexFront(int x, int y, int z);

   /**
    * Removes all vertices from conformation. After calling this function,
    * vertices should not be accessed until new vertices are added.
    */
   virtual void clear();

   /**
    * Writes conformation as a single space separated list of coordinates.
    * @param os
    */
   virtual void writeAsText(std::ostream& os) const;

   /**
    * Gets conformation as a string in same formate as that output by
    * writeAsText(std::ostream& os) const.
    * @return a string containing coordinates separated by single spaces.
    */
   virtual std::string writeAsText() const;

   // TODO:
   //  virtual void writeAsUofS(std::ostream& os) const;
   //  virtual void writeAsNewsud(std::ostream& os) const;

   // TODO: should throw exception for io error
   /**
    * Outputs conformation in Rob Scharein's cube format. This formation is a
    * binary representation of news format where each vertex is written in three
    * bits.
    * @param os
    * @return false for io exception
    */
   virtual bool writeAsCube(std::ostream& os) const;

   /**
    * Reads one conformation written as a line in a text file where coordinates
    * are separated by spaces. If line contains 3*n+1 numbers, assumes first
    * number is number of vertices and ignores that number.  
    * @param is
    * @return false for io exception or end of stream.
    */
   virtual bool readFromText(std::istream& is);

   /**
    * Reads one conformation from a coord file, that is a text file where
    * each vertex is written as three space separated integers, one vertex per
    * line. The end of conformation is denoted by a line with less than three
    * integers.  
    * @param is
    * @return false for io exception or end of stream.
    */
   virtual bool readFromCoords(std::istream& is);

   /**
    * Parses a string in same format as one line read by 
    * readFromText(std::istream& is).
    * @param text assumed to be a string of representation of integers separated
    *  by spaces.
    * @return 
    */
   virtual bool readFromText(const std::string& text);

   // TODO:
   //  virtual bool readFromUofS(std::istream& is);
   //  virtual bool readFromNewsud(std::istream& is);

   /**
    * Inputs conformation in Rob Scharein's cube format. This formation is a
    * binary representation of news format where each vertex is written in three
    * bits.
    * @param os
    * @return false for io exception
    */
   virtual bool readFromCube(std::istream& is);

   /**
    * Translates every vertex of conformation by same offset.
    * @param offset added to every vertex of conformation.
    */
   virtual void translate(const threevector<int>& offset);

   /**
    * Translates every vertex of conformation by same offset.
    * @param x added to x-coordinate of every vertex of conformation.
    * @param y added to y-coordinate of every vertex of conformation.
    * @param z added to z-coordinate of every vertex of conformation.
    */
   virtual void translate(int x, int y, int z);

   /**
    * Computes radius of gyration of conformation.
    * @return radius of gyration.
    */
   virtual double rog() const;

   /**
    * Computes writhe and ACN.
    * @param writhe will be set to writhe of conformation.
    * @param acn will be set to ACN of conformation.
    */
   virtual void writheACN(double& writhe, double& acn) const;

   /**
    * Computes writhe and ACN between two conformations.
    * @param other should not intersect this at any edge. 
    * @param writhe will be set to writhe of conformation.
    * @param acn will be set to ACN of conformation.
    */
   virtual void writheACN(const clkConformationAsList& other, double& writhe, double& acn) const;


protected:

   // TODO: make template class.
   std::deque<threevector<int> > data;

private:
   bool lineIsVertex(std::string line);

};

#endif	/* CLKCONFORMATIONASLIST_H */

