#pragma once

#include "clkTypes.h"

bool add_edge_to_knot(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end);
void bfacf_init_pool_and_lattice(CubicLatticeKnotPtr clkp, bool ignore_mid, int poolsize, ivector min, ivector max);
void clk_get_statistics(CubicLatticeKnotPtr clkp, int& n1, int& n2, int& n3, int& n4);
bool clkp_always_turns(CubicLatticeKnotPtr clkp);
