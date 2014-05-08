/* 
 * File:   bfacfProbabilitiesFromZ.cpp
 * Author: kmo
 * 
 * Created on January 19, 2013, 3:19 PM
 */

#include "bfacfProbabilitiesFromZ.h"

bfacfProbabilitiesFromZ::bfacfProbabilitiesFromZ() { }

bfacfProbabilitiesFromZ::bfacfProbabilitiesFromZ(const bfacfProbabilitiesFromZ& orig) { }

bfacfProbabilitiesFromZ::~bfacfProbabilitiesFromZ() { }

// TODO: Find correct source for probabilities.

double bfacfProbabilitiesFromZ::p2() const
{
   double zSquared = z() * z();
   return zSquared * m2();
}

double bfacfProbabilitiesFromZ::p0() const
{
//   return 0.5 - p_plus2();
   return 0.5* (p2()+m2());
}

double bfacfProbabilitiesFromZ::m2() const
{
   double zSquared = z() * z();
//   return zSquared / (1.0 + 3.0 * zSquared);
   return 1.0 / (1.0 + 3.0 * zSquared);
}

double bfacfProbabilitiesFromZ::p4p2() const
{
   return 4.0 * p2();
}

double bfacfProbabilitiesFromZ::p03p2() const
{
   return p0() + 3.0 * p2();
}

double bfacfProbabilitiesFromZ::m23p2() const
{
   return m2() + 3.0 * p2();
}

double bfacfProbabilitiesFromZ::p2p0() const
{
   return 2.0 * p0();
}

double bfacfProbabilitiesFromZ::p2p02p2() const
{
   return 2.0 * p0() + 2.0 * p2();
}

