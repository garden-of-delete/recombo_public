#pragma once

#include "clkTypes.h"

void swapEdgesPool(CubicLatticeKnotPtr clkp, int loc1, int loc2);
void freezeEdge(CubicLatticeKnotPtr clkp, EdgePtr ep);
bool freezeEdge(CubicLatticeKnotPtr clkp, ivector start);
void freezeEdge(EdgePtr ep);
void thawEdge(CubicLatticeKnotPtr clkp, EdgePtr ep);
void thawEdge(EdgePtr ep);
bool freeze_component(CubicLatticeKnotPtr clkp, int ID, bool freeze);

extern int clk_minedges, clk_maxedges;
void clk_set_edge_limit(int min, int max);
void clk_set_nedge_limit(CubicLatticeKnotPtr clkp, int min, int max, int ID);
