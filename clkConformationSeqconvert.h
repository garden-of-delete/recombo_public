/* 
 * File:   clkConformationSeqconvert.h
 * Author: kmo
 *
 * Created on December 4, 2012, 3:08 PM
 */

#ifndef CLKCONFORMATIONSEQCONVERT_H
#define	CLKCONFORMATIONSEQCONVERT_H

#include "clk.h"

class clkConformationSeqconvert : public clk
{
public:
   clkConformationSeqconvert();
   clkConformationSeqconvert(const clk& orig);
   clkConformationSeqconvert(const clkConformationSeqconvert& orig);
   virtual ~clkConformationSeqconvert();

   virtual int size() const;
   virtual void getVertex(int index, threevector<int>& v) const;

   //  virtual void putKnot(FILE *fp) = 0;
   //  virtual bool getKnot(FILE *fp) = 0;

   virtual void computeWritheACN(double& writhe, double& ACN);
   virtual double radiusOfGyration();

   virtual clkConformationSeqconvert& operator=(const clk& orig);
   virtual clkConformationSeqconvert& operator=(const clkConformationSeqconvert& orig);

private:

   class impl;
   impl* implementation;
};

#endif	/* CLKCONFORMATIONSEQCONVERT_H */

