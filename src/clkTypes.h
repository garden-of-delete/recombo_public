#pragma once

#include "legacy.h"

typedef struct _edge
{
    struct _edge* prev; // pointer to previous edge in doubly linked list of edges
    struct _edge* next; // pointer to next edge in doubly linked list of edges
    int dir; // direction of edge, integer from 0 to 5, inclusive
    int increment[3]; // increment corresponding to dir above
    int start[3]; // position of start of edge, a 3-tuple of integers
    int locpool; // location of this edge in edge pool
    bool frozen; // whether is frozen
    void* comp; // component this edge is located in
    int colour; // edge colour
    int ID; // unique ID (for debugging purposes), assigned when edge created
    int scratch; // temporary space
} Edge;

typedef Edge* EdgePtr;

typedef struct
{
    EdgePtr ep1, ep2;
} EdgePair;

typedef EdgePair* EdgePairPtr;

#define COMPONENT_CLK_FLAG_OPEN  1  // open-ended

typedef struct _comp
{
    int ID;
    int flags; // properties of this component
    int nedges; // number of edges in knot component
    int minedges; // minimum number of edges permitted in this component
    int maxedges; // maximum number of edges peymitted in this component
    void* clkp; // which clkp this component is in
    EdgePtr first_edge; // pointer to first edge in knot component
    EdgePtr last_edge; // pointer to last edge in knot component
    struct _comp* prev; // previous component
    struct _comp* next; // next component

    double z; // fugacity parameter, might be different for different knots

    // values of the next three parameters computed from z value
    double p_minus2; // p(-2), probability of doing a -2 move
    double p_0; // p(0), probability of doing a 0 move
    double p_plus2; // p(+2), probability of doing a +2 move

    // the following are precomputed in hope of an improvement in running time
    double p_4p2; // 4 * p(+2)
    double p_03p2; // p(0) + 3 * p(+2)
    double p_m23p2; // p(-2) + 3 * p(+2)
    double p_2p0; // 2 * p(0)
    double p_2p02p2; // 2 * p(0) + 2 * p(+2)
} ComponentCLK; // cubic lattice knot component

typedef ComponentCLK* ComponentCLKPtr;

typedef struct
{
    int loffset[3]; // offset of first vertex in self-avoidance lattice
    int nedges_total; // total number of edges in knot (note: not number in edge pool)
    int nedges_alt; // not used except by KnotPlot (for now)

    EdgePtr* edgepool; // pool of pointers to edges (array of size `poolsize')
    int poolsize; // maximum number of edges in pool
    int nfrozen; // number of `frozen edges', frozen edges are never chosen
    // (this allows for open-ended strings or the inclusion
    // of topological obstructions)
    bool auto_recentre; // whether or not to recentre after an edge hit

    char* lattice; // lattice used for self-avoidance checking
    int* alt_lattice; // alternate lattice, allocated on demand

    ivector max_range; // maximum range allowed

    ComponentCLKPtr fcomp; // first cubic lattice knot component in a linked list
    ComponentCLKPtr lcomp; // last cubic lattice knot component in a linked list
    int ncomps; // number of components


    // the following are to keep track of various things that may or may not be of interest
    // but are not actually used in the algorithm
    int edge_hits[6]; // number of times knot hits edge of self-avoidance lattice
    int success_minus2; // number of successful -2 moves
    int success_0; // number of successful 0 moves
    int success_plus2; // number of successful +2 moves
    int p2dir[6];
    int pm2dir[6];
    int total_sampled;
    int total_output;
    int wrong_length;

    FILE* fpout;

    int iteration; // iteration count
} CubicLatticeKnot;

typedef CubicLatticeKnot* CubicLatticeKnotPtr;
