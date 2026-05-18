#include "clkInit.h"
#include "clkLattice.h"
#include "clkConstants.h"
#include "clkTables.h"
#include "clkGeometry.h"
#include "clkValidation.h"
#include "bfacfMove.h"
#include "edgePool.h"

#include <cstdio>
#include <cstdlib>

bool add_edge_to_knot(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end)
{
    EdgePtr ep = (EdgePtr)calloc(1, sizeof(Edge)); // create a new edge

    sub_ivector(ep->increment, end, start); // increment is (end - start)
    ep->dir = clk_direction(ep->increment); // direction associated with increment
    copy_ivector(ep->start, start); // starting location of this edge
    ep->frozen = false; // edge is free to move
    ep->comp = (void*)comp;

    if (comp->first_edge)
    {
        // some edges have already been added
        comp->last_edge->next = ep;
        ep->prev = comp->last_edge;
        comp->last_edge = ep;
    }
    else
    {
        // this is the first edge
        comp->first_edge = comp->last_edge = ep;
    }

    comp->first_edge->prev = comp->last_edge;
    comp->last_edge->next = comp->first_edge;
    comp->nedges++;
    clkp->nedges_total++;

    return clk_check_increment(ep->increment);
}

void bfacf_init_pool_and_lattice(CubicLatticeKnotPtr clkp, bool ignore_mid, int poolsize, ivector min, ivector max)
{
    ComponentCLKPtr comp = clkp->fcomp;

    clkp->auto_recentre = true;

    // create a `swimming space' (a pool) for our edges
    // this is to avoid frequent memory allocation / deallocation

    clkp->poolsize = poolsize;
    clkp->edgepool = (EdgePtr*)calloc(clkp->poolsize, sizeof(EdgePtr));
    int index = 0;
    int ID = 0;

    comp = clkp->fcomp;
    while (comp)
    {
        comp->ID = ID++;
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (index > poolsize)
            {
                printf("yikes! %d %d\n", index, poolsize);
                exit(11);
            }
            clkp->edgepool[index++] = ep;
            ep = ep->next;
        }
        comp = comp->next;
    }

    while (index < poolsize)
        clkp->edgepool[index++] = (EdgePtr)calloc(1, sizeof(Edge)); // new, unused edge

    for (index = 0; index < clkp->poolsize; index++)
    {
        clkp->edgepool[index]->locpool = index;
        clkp->edgepool[index]->ID = index;
    }

    clkp->lattice = (char*)calloc(LATTICE_TOTAL_SIZE, sizeof(char));
    ivector mid;
    midpoint(mid, min, max);

    if (ignore_mid)
        set_ivector(mid, 0, 0, 0);

    init_lattice(clkp, mid);


    comp = clkp->fcomp;
    while (comp)
    {
        ivector closing_incr;
        sub_ivector(closing_incr, comp->first_edge->start, comp->last_edge->start);
        // if the edge that would close the knot is not valid,
        // then we have an open-ended string.
        if (!clk_check_increment(closing_incr))
        {
            comp->flags |= COMPONENT_CLK_FLAG_OPEN;
            delete_Edge(clkp, comp, comp->last_edge);
            freezeEdge(clkp, comp->first_edge);
            freezeEdge(clkp, comp->last_edge);
            comp->first_edge->frozen = comp->last_edge->frozen = true;
        }

        comp = comp->next;
    }
}


void clk_get_statistics(CubicLatticeKnotPtr clkp, int& n1, int& n2, int& n3, int& n4)
{
    // gather statistics on local edge configurations with an eye to possibly
    // reordering the cases
    // currently, this function is used only by the KnotPlot driver
    if (!clkp) return;
    n1 = n2 = n3 = n4 = 0;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            // Case 1: chosen edge is parallel to both adjacent edges,
            // all four moves are +2 moves.

            // directions of adjacent edges same as chosen edge

            if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
            {
                ++n1;
            }

            // Case 2: chosen edge is perpendicular to both adjacent edges
            // and those adjacent edges are anti-parallel to each other.
            // One move is a -2 move, other three are +2 moves.

            else if (perp[ep->dir][ep->next->dir] && anti[ep->prev->dir][ep->next->dir])
            {
                ++n2;
            }

            // Case 3: chosen edge is perpendicular to one adjacent edge
            // and parallel to the other adjacent edge.
            // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
            // because the case of both being parallel has been ruled out above.
            // One move is a 0 move, other three are +2 moves.

            // NOTE: test this before Case 1 above

            else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
            {
                ++n3;
            }

            // Case 4: Only possibility remaining: chosen edge is perpendicular
            // to both adjacent edges, and the adjacent edges are not anti-parallel.
            // Two moves are 0 moves and two moves are +2 moves.

            else
            {
                ++n4;
            }

            ep = ep->next;
        }
        comp = comp->next;
    }
}

bool clkp_always_turns(CubicLatticeKnotPtr clkp)
{
    if (!clkp) return false;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (ep->dir == ep->next->dir) return false;
            ep = ep->next;
        }
        comp = comp->next;
    }
    return true;
}
