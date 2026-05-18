#pragma once

#include "clkTypes.h"

// Lookup tables for BFACF moves
extern int opposite[6];
extern int increment_NEWSUD[6][3];
extern int turn[6][4];
extern int reflect[6][6];
extern int kross[6][6];
extern int rotate90[6][6];
extern bool perp[6][6];
extern bool anti[6][6];
extern char clk_dir_name[6][6];

// Move functions for BFACF algorithm
bool perform_move(CubicLatticeKnotPtr knot);

// Helper functions for specific move types
bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);
bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep, ivector increment, int dir);
bool perform_plus2_move_alt(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);
bool perform_0_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);
bool perform_minus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);

// Helper function
void delete_Edge(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);
