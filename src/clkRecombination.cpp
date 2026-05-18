#include "clkRecombination.h"
#include "bfacf.h"
#include "clkConstants.h"
#include "clkTables.h"
#include "clkGeometry.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>

extern int clk_count_edges(EdgePtr ep);
extern void pexit(char* s);

bool clk_recombo_limit = false;
int clk_min_arc_length_distance = 0;
bool clk_recombo_direct = true;

void clk_fix_incr(EdgePtr ep)
{
    sub_ivector(ep->increment, ep->next->start, ep->start);
    ep->dir = clk_direction(ep->increment);
    if (ep->dir == MOVE_INVALID)
        pexit("invalid increment detected");
}

bool perform_recombination_inverted(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
    if (ep1->comp != ep2->comp) return false;
    ivector test;
    sub_ivector(test, ep1->start, ep2->start);

    // the following section can be removed once we use pairs
    if (!clk_check_increment(test))
        pexit("impossible situation perform_recombination_inverted () A");
    if (ep1->next == ep2 || ep1->prev == ep2)
        pexit("impossible situation perform_recombination_inverted () B");

    EdgePtr ep1n = ep1->next;
    EdgePtr ep2n = ep2->next;
    EdgePtr ep1nn = ep1n->next;
    EdgePtr ep2p = ep2->prev;

    EdgePtr epp = ep1n;
    EdgePtr ep = ep1n->next;

    int kount = 0;

    while (ep != ep2)
    {
        EdgePtr epn = ep->next;
        ep->next = epp;
        ep->prev = epn;
        epp = ep;
        ep = epn;
        if (++kount > 16008)
            pexit("infinite loop A");
    }

    ep1->next = ep2;
    ep2->prev = ep1;
    ep2->next = ep2p;
    ep1n->prev = ep1nn;
    ep1n->next = ep2n;
    ep2n->prev = ep1n;

    kount = 0;
    ep = ep1;
    do
    {
        clk_fix_incr(ep);
        ep = ep->next;
        if (++kount > 16008)
            pexit("infinite loop B");
    }
    while (ep != ep2n);

    clk_validate(clkp, "end of perform_recombination_inverted()");

    return true;
}

bool perform_recombination(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
    if (!ep1 || !ep2)
    {
        fprintf(stderr, " *** perform_recombination (): null pointer");
        if (!ep1) fprintf(stderr, "ep1 is null ");
        if (!ep2) fprintf(stderr, "ep2 is null ");
        fflush(stderr);

        return false;
    }
    // following is needed to prevent edges like 29 and 31 in 9jun10a.k from being recombined
    if (ep1->next == ep2->prev || ep1->prev == ep2->next) return false;

    if (ep1->dir == ep2->dir)
        return perform_recombination_inverted(clkp, ep1, ep2);

    // do some paranoia checking (for now)
    if (!anti[ep1->dir][ep2->dir])
        pexit("edges are not anti-parallel as expected!");


    ivector incr;
    sub_ivector(incr, ep1->next->start, ep2->start);
    if (!clk_check_increment(incr))
        pexit("invalid increment found!");

    ComponentCLKPtr comp1 = (ComponentCLKPtr)ep1->comp;
    ComponentCLKPtr comp2 = (ComponentCLKPtr)ep2->comp;

    if (comp1 == comp2 && comp1->nedges < 5) return false;
    if (clk_recombo_limit && comp1 == comp2 && clkp->ncomps > 1) return false;

    EdgePtr ep1n = ep1->next;
    EdgePtr ep1p = ep1->prev;
    EdgePtr ep2n = ep2->next;
    EdgePtr ep2p = ep2->prev;

    // update directions of edges
    sub_ivector(ep1->increment, ep2->next->start, ep1->start);
    ep1->dir = clk_direction(ep1->increment);
    sub_ivector(ep2->increment, ep1->next->start, ep2->start);
    ep2->dir = clk_direction(ep2->increment);

    // update how things are connected
    ep1->next = ep2n;
    ep2->next = ep1n;
    ep1n->prev = ep2;
    ep2n->prev = ep1;

    if (comp1 == comp2)
    {
        comp2 = (ComponentCLKPtr)calloc(1, sizeof(ComponentCLK));
        comp2->minedges = clk_minedges;
        comp2->maxedges = clk_maxedges;
        memcpy(comp2, comp1, sizeof(ComponentCLK));
        comp1->first_edge = ep1;
        comp1->last_edge = ep1->prev;
        comp2->first_edge = ep2;
        comp2->last_edge = ep2->prev;
        comp1->nedges = clk_count_edges(ep1);
        comp2->nedges = clk_count_edges(ep2);
#ifdef MACOSX
        if (comp1->nedges < 4 || comp2->nedges < 4) system("say oh dear, not again &");
#endif
        // add new component comp2 to end of linked list
        comp2->prev = clkp->lcomp;
        comp2->next = (ComponentCLKPtr)NULL;
        clkp->lcomp->next = comp2;
        clkp->lcomp = comp2;
        clkp->ncomps++;
    }
    else
    {
        comp1->first_edge = ep1;
        comp1->last_edge = ep1->prev;
        comp1->nedges = clk_count_edges(ep1);

        // remove comp2 from linked list
        if (comp2->prev)
            comp2->prev->next = comp2->next;
        else
            clkp->fcomp = comp2->next;

        if (comp2->next)
            comp2->next->prev = comp2->prev;
        else
            clkp->lcomp = comp2->prev;

        free(comp2);

        clkp->ncomps--;
    }

    // recompute IDs
    ComponentCLKPtr comp = clkp->fcomp;
    int ID = 0;
    while (comp)
    {
        comp->ID = ID++;
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            ep->comp = (void*)comp;
            ep = ep->next;
        }
        comp = comp->next;
    }

    //  validate (clkp, "end of perform_recombination()");

    return true;
}

bool clk_allocation_alt_lattice(CubicLatticeKnotPtr clkp)
{
    if (!clkp) return false;
    if (!clkp->alt_lattice)
        clkp->alt_lattice = (int*)calloc(LATTICE_TOTAL_SIZE, sizeof(int));
    return true;
}

int bfacf_perform_recombination(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr))
{
    if (!clk_allocation_alt_lattice(clkp)) return 0;


    //  validate (clkp, "start of A");

    // fill knot path with edge pool location
    ComponentCLKPtr comp = clkp->fcomp;
    EdgePtr ep1;
    while (comp)
    {
        if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
        {
            ep1 = comp->first_edge;
            for (int i = 0; i < comp->nedges; i++)
            {
                clkp->alt_lattice[lat(ep1->start)] = ep1->locpool + 1;
                ep1->frozen = false;
                ep1 = ep1->next;
            }
        }
        comp = comp->next;
    }

    // now go thru the edges looking for antiparallel edges at distance 1
    comp = clkp->fcomp;
    int kount = 0;
    // bfacf_lattice_info (clkp);

    while (comp)
    {
        if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
        {
            ep1 = comp->first_edge;
            for (int i = 0; i < comp->nedges; i++)
            {
                for (int dir = 0; dir < 4; dir++)
                {
                    ivector test_loc;
                    add_ivector(test_loc, ep1->start, increment_NEWSUD[turn[ep1->dir][dir]]);
                    int IDp1 = clkp->alt_lattice[lat(test_loc)];
                    if (IDp1)
                    {
                        if (IDp1 < 1 || IDp1 >= clkp->poolsize)
                        {
                            printf("bad value: %d\n", IDp1);
                            exit(22);
                        }
                        EdgePtr ep2;
                        bool ok = true;
                        if (clk_recombo_direct)
                        {
                            ep2 = clkp->edgepool[IDp1 - 1]->prev;
                            if (clk_recombo_limit && ep1->comp == ep2->comp && clkp->ncomps > 1)
                                ok = false;
                            if (ok && anti[ep1->dir][ep2->dir] && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2)
                                >= clk_min_arc_length_distance &&
                                !(ep1->next == ep2->prev || ep1->prev == ep2->next))
                            {
                                ++kount;
                                mark(ep1, ep2);
                            }
                        }
                        else
                        {
                            ep2 = clkp->edgepool[IDp1 - 1];
                            if (ep1->dir == ep2->dir && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2) >=
                                clk_min_arc_length_distance &&
                                !(ep1->next == ep2->prev || ep1->prev == ep2->next))
                            {
                                ++kount;
                                mark(ep1, ep2);
                            }
                        }
                    }
                }
                ep1 = ep1->next;
            }
        }
        comp = comp->next;
    }

    comp = clkp->fcomp;
    while (comp)
    {
        if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
        {
            ep1 = comp->first_edge;
            for (int i = 0; i < comp->nedges; i++)
            {
                clkp->alt_lattice[lat(ep1->start)] = 0;
                ep1 = ep1->next;
            }
        }
        comp = comp->next;
    }

    return kount;
}
