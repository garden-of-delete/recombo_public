#pragma once

#include "clkTypes.h"

extern int opposite[6];
extern ivector increment_NEWSUD[6];
extern int turn[6][4];
extern bool perp[6][6];
extern bool anti[6][6];
extern char clk_dir_name[6][6];
extern int kross[6][6];

extern int clk_minedges, clk_maxedges;
extern bool clk_recombo_limit, clk_recombo_direct;
extern int clk_min_arc_length_distance;
