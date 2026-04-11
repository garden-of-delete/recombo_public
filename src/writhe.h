#ifndef WRITHE_H
#define WRITHE_H

#include "legacy.h"

double writhe(double &ACN, int nvert, ivector *vert, double jitter);
double writhe_open(ivector *AB, ivector *CD);
double radius_of_gyration(int nvert, ivector *vert);

#endif
