/* 
 * File:   clkCigar.h
 * Author: kmo
 *
 * Created on December 7, 2012, 7:47 PM
 */

#ifndef CLKCIGAR_H
#define	CLKCIGAR_H

#include "clk.h"

/**
 * Implements an unknotted conformation of fixed length in the cubic lattice.
 * The geometry of this kind of conformation will depend only on length.
 */
class clkCigar : public clk
{
public:
   /**
    * Default constructor produces a square whose first vertex is the origin.
    */
   clkCigar();
   
   /**
    * Constructs a
    * @param length the number of vertices desired. Must be an even integer 
    * greater than four.
    */
   clkCigar(int length);
   
   /**
    * Copy constructor. Will produce an unknotted conformation with same number
    * of vertices.
    * @param orig
    */
   clkCigar(const clkCigar& orig);
   
   /**
    * Destructor.
    */
   virtual ~clkCigar();

   virtual int size() const;
   virtual void getVertex(int index, threevector<int>& v) const;

private:

   const int halfLength;
};

#endif	/* CLKCIGAR_H */

