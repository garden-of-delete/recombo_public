/* 
 * File:   bfacfProbabilitiesFromZ.h
 * Author: kmo
 *
 * Created on January 19, 2013, 3:19 PM
 */

#ifndef BFACFPROBABILITIESFROMZ_H
#define	BFACFPROBABILITIESFROMZ_H

#include "bfacfProbabilities.h"

/**
 * Class specifying bfacf probabilities calculated from z-value.
 */
class bfacfProbabilitiesFromZ: public bfacfProbabilities
{
public:
   bfacfProbabilitiesFromZ();
   bfacfProbabilitiesFromZ(const bfacfProbabilitiesFromZ& orig);
   virtual ~bfacfProbabilitiesFromZ();

   virtual double p2() const;
   virtual double p0() const;
   virtual double m2() const;
   virtual double p4p2() const;
   virtual double p03p2() const;
   virtual double m23p2() const;
   virtual double p2p0() const;
   virtual double p2p02p2() const;
   
   /**
    * Accessor for z-value.
    * @return The z-value.
    */
   virtual double z() const = 0;
   
private:

};

#endif	/* BFACFPROBABILITIESFROMZ_H */

