/* 
 * File:   clkCigar.cpp
 * Author: kmo
 * 
 * Created on December 7, 2012, 7:47 PM
 */

#include "clkCigar.h"
#include "threevector.h"

clkCigar::clkCigar(): halfLength(2) { }

clkCigar::clkCigar(int length): halfLength(length/2) { }

clkCigar::clkCigar(const clkCigar& orig): clk(orig), halfLength(orig.halfLength)
{ }

clkCigar::~clkCigar() { }

int clkCigar::size() const
{
   return 2 * halfLength;
}


void clkCigar::getVertex(int index, threevector<int>& v) const
{
   int x, y;
   if((index / halfLength) & 1)
   {
      // index / halfLength is odd
      x = halfLength - 1 - (index % halfLength);
      y = 1;
   }
   else
   {
      // index / halfLength is even
      x = index % halfLength;
      y = 0;
   }
   
   v.set(x, y, 0);
}
