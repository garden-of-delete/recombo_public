/* 
 * File:   regularNgon.h
 * Author: kmo
 *
 * Created on November 30, 2012, 6:49 PM
 */

#ifndef REGULARNGON_H
#define	REGULARNGON_H

#include "conformation.h"

class regularNgon: public conformation
{
public:
   regularNgon();
   regularNgon(int numberSides, double radius);
   regularNgon(const regularNgon& orig);
   virtual ~regularNgon();
   
   virtual int size() const;
   virtual void getVertex(int index, threevector<double>& v) const;
      
   virtual double getRadius() const;
   
protected:
   const int numberSides;
   const double radius;
};

#endif	/* regularNgon_H */