/* 
 * File:   regularNgon.cpp
 * Author: kmo
 * 
 * Created on November 30, 2012, 6:49 PM
 */

#include "regularNgon.h"
#define M_PI 3.14159265358979323846 //added for windows compatibility
#include <cmath>

using namespace std;

regularNgon::regularNgon() : conformation(), numberSides(3), radius(1) { }

regularNgon::regularNgon(int numberSides, double radius) : conformation(),
numberSides(numberSides), radius(radius) { }

regularNgon::regularNgon(const regularNgon& orig) : conformation(orig),
numberSides(orig.numberSides), radius(orig.numberSides) { }

regularNgon::~regularNgon() { }

int regularNgon::size() const
{
   return numberSides;
}

void regularNgon::getVertex(int index, threevector<double>& v) const
{
   v.set(radius * cos((2.0 * M_PI * index) / numberSides), radius * sin((2.0 * M_PI * index) / numberSides), 0.0);
}

double regularNgon::getRadius() const
{
   return radius;
}

