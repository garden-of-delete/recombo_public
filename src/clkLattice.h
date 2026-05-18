#pragma once

#include "clkTypes.h"

void init_lattice(CubicLatticeKnotPtr clkp, ivector mid);
bool recentre_knot_in_lattice(CubicLatticeKnotPtr knot);

void set_lattice(CubicLatticeKnotPtr clkp, char value, ivector min, ivector max);
void set_lattice(CubicLatticeKnotPtr clkp, char value);

void fill_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max);
void fill_lattice(CubicLatticeKnotPtr clkp);

void clear_lattice(CubicLatticeKnotPtr clkp);
void clear_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max);

void invert_lattice(CubicLatticeKnotPtr clkp);
void dilate_lattice(CubicLatticeKnotPtr clkp);
int bfacf_info_filled(CubicLatticeKnotPtr knot);
