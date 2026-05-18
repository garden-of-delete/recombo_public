#pragma once

#include "clkTypes.h"

void clk_fix_incr(EdgePtr ep);
bool perform_recombination_inverted(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2);
bool perform_recombination(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2);
bool clk_allocation_alt_lattice(CubicLatticeKnotPtr clkp);
int bfacf_perform_recombination(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr));

extern bool clk_recombo_limit, clk_recombo_direct;
extern int clk_min_arc_length_distance;
