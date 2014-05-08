/* 
 * File:   clk.h
 * Author: kmo
 *
 * Created on November 27, 2012, 6:49 PM
 */

#ifndef CLK_H
#define	CLK_H

#include <ostream>

#include "conformation.h"

/**
 * Class representing conformations whose vertices have integer coordinates.
 */
class clk : public conformation
{
public:
   /**
    * Default constructor.
    */
   clk();
   
   /**
    * Copy constructor.
    * @param orig
    */
   clk(const clk& orig);
   
   /**
    * Destructor.
    */
   virtual ~clk();

   /**
    * Populates a list with coordinates of vertices.
    * @param list Will be cleared before being populated with vertex coordinates.
    */
   virtual void getVertex(int n, threevector<int>& v) const = 0;
   
   virtual void getVertex(int n, threevector<double>& v) const;
   
    /**
    * Populates a list with coordinates of vertices.
    * @param list Will be cleared before being populated with vertex coordinates.
    */
   virtual void getVertices(std::list<threevector<int> >&) const;
   
   /**
    * Tests two conformations for equality.
    * @param other conformation.
    * @return false if *this and other have different length or differ at
    * any vertex, and true otherwise.
    */
   virtual bool operator==(const clk& other) const;

   /**
    * Tests two conformations for inequality.
    * @param other conformation.
    * @return true if *this and other have different length or differ at
    * any vertex, and false otherwise.
    */
   virtual bool operator!=(const clk& other) const;

protected:
   virtual void display(std::ostream& out) const;

};

#endif	/* CLK_H */

