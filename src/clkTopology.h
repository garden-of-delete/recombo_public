#pragma once

#include "clkTypes.h"

// Topology detection
EdgePtr is_tight_clasp(CubicLatticeKnotPtr clkp, EdgePtr ep);
EdgePtr is_parsite(CubicLatticeKnotPtr clkp, EdgePtr ep);
bool has_tight_clasp(CubicLatticeKnotPtr clkp);
bool has_parsite(CubicLatticeKnotPtr clkp);
int numb_tight_clasps(CubicLatticeKnotPtr clkp);
int numb_tight_clasps(CubicLatticeKnotPtr clkp, void mark(ivector, ivector));
int numb_parsites(CubicLatticeKnotPtr clkp);
int numb_parsites(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr));

// Topology moves
bool bfacf_strand_pass(CubicLatticeKnotPtr clkp);
bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, EdgePtr ep);
bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp);
