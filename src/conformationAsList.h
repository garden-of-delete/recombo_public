/* 
 * File:   conformationAsList.h
 * Author: kmo
 *
 * Created on December 4, 2012, 3:04 PM
 */

#ifndef CONFORMATIONASLIST_H
#define	CONFORMATIONASLIST_H

#include <list>
#include <vector>

#include "conformation.h"

class conformationAsList : public conformation
{
public:
   conformationAsList();
   conformationAsList(const conformationAsList& orig);
   conformationAsList(const conformation& orig);
   virtual ~conformationAsList();

   virtual int size() const;
   virtual void getVertex(int index, threevector<double>& v) const;

   /**
    * Add a new vertex in new last position of conformation. This operation will
    * increase size by one.
    * @param v vertex to add.
    */
   virtual void addVertexBack(const threevector<double>& v);

   /**
    * Add a new vertex in new last position of conformation. This operation will
    * increase size by one.
    * @param x first coordinate of vertex to add.
    * @param y second coordinate of vertex to add.
    * @param z third coordinate of vertex to add.
    */
   virtual void addVertexBack(double x, double y, double z);

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

   /**
    * Reads one conformation written as a line in a text file where coordinates
    * are separated by spaces. If line contains 3*n+1 numbers, assumes first
    * number is number of vertices and ignores that number.  
    * @param is
    * @return false for io exception or end of stream.
    */
   virtual bool readFromText(std::istream& is);

   /**
    * Parses a string in same format as one line read by 
    * readFromText(std::istream& is).
    * @param text assumed to be a string of representation of integers separated
    *  by spaces.
    * @return 
    */
   virtual bool readFromText(const std::string& text);

   /**
    * Translates every vertex of conformation by same offset.
    * @param offset added to every vertex of conformation.
    */
   virtual void translate(const threevector<double>& offset);

   /**
    * Translates every vertex of conformation by same offset.
    * @param x added to x-coordinate of every vertex of conformation.
    * @param y added to y-coordinate of every vertex of conformation.
    * @param z added to z-coordinate of every vertex of conformation.
    */
   virtual void translate(double x, double y, double z);

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

protected:
   // TODO: make template class.
   std::vector<threevector<double> > data;

};

#endif	/* CONFORMATIONASLIST_H */

