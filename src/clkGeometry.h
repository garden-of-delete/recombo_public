#pragma once

#include "clkTypes.h"

EdgePtr clk_get_edge(CubicLatticeKnotPtr clkp, ivector start);
int clk_count_edges(EdgePtr ep);
int clk_arc_length_distance(EdgePtr ep1, EdgePtr ep2);
int clk_arc_length_distance(CubicLatticeKnotPtr clkp, ivector s1, ivector s2);
void clk_get_extent(ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot);
void clk_bounding_box(int& X, int& Y, int& Z, CubicLatticeKnotPtr knot);
int clk_direction(ivector incr);
int clk_direction(ivector end, ivector start);
