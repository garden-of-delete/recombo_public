/* 
 * File:   bfacfProbabilitiesFromZFixed.h
 * Author: kmo
 *
 * Created on January 19, 2013, 3:19 PM
 */

#ifndef BFACFPROBABILITIESFROMZFIXED_H
#define	BFACFPROBABILITIESFROMZFIXED_H

#include "bfacfProbabilitiesFromZ.h"

/**
 * Probabilities for bfacf moves, computed from z-value which is fixed at time
 * of construction.
 */
class bfacfProbabilitiesFromZFixed: public bfacfProbabilitiesFromZ
{
public:
   bfacfProbabilitiesFromZFixed();
   bfacfProbabilitiesFromZFixed(double z);
   bfacfProbabilitiesFromZFixed(const bfacfProbabilitiesFromZ& orig);
   virtual ~bfacfProbabilitiesFromZFixed();

   virtual double z() const;
   
private:
   const double zValue;
};

#endif	/* BFACFPROBABILITIESFROMZFIXED_H */

