#pragma once

#include "clkTypes.h"

// Structure validation
int clk_validate(CubicLatticeKnotPtr clkp);
void clk_validate(CubicLatticeKnotPtr clkp, char* s);

// Increment validation
bool clk_check_increment(ivector incr);
void clk_check_increments(CubicLatticeKnotPtr knot);

// Lattice validation
bool clk_check_lattice(CubicLatticeKnotPtr knot);
bool clk_check_for_edge_hits(CubicLatticeKnotPtr knot, int dir, ivector test_location);
