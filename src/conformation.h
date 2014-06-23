/* 
 * File:   conformation.h
 * Author: kmo
 *
 * Created on November 30, 2012, 6:49 PM
 */

#ifndef CONFORMATION_H
#define	CONFORMATION_H

#include <ostream>
#include <list>

#include "threevector.h"

/**
 * Class representing polygonal paths.
 */
class conformation
{
public:
   /**
    * Default constructor.
    */
   conformation();

   /**
    * Copy constructor.
    */
   conformation(const conformation& orig);

   /**
    * Destructor.
    */
   virtual ~conformation();

   /**
    * Accesses length of conformation.
    * @return The number of vertices in conformation.
    */
   virtual int size() const = 0;
   
   /**
    * Accesses coordinates of one vertex of conformation.
    * @param n Position of vertex to access.
    * @param v Coordinates of v will be set to the coordinates of nth vertex.
    */
   virtual void getVertex(int n, threevector<double>& v) const = 0;
   
   /**
    * Populates a list with coordinates of vertices.
    * @param list Will be cleared before being populated with vertex coordinates.
    */
   virtual void getVertices(std::list<threevector<double> >& list) const;

   /**
    * Displays coordinates on a single line, each coordinate separated by a 
    * single space. Does not append newline or any other character.
    * @param out an ostream.
    */
   friend std::ostream& operator<<(std::ostream& out, const conformation& conformation);

protected:
   
   /**
    * Function called by operator<<(std::ostream&, const conformation&).
    * @param out
    */
   virtual void display(std::ostream& out) const;
};

#endif	/* CONFORMATION_H */

