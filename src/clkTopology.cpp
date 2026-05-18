#include "clkTopology.h"
#include "clkConstants.h"
#include "clkTables.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

EdgePtr is_tight_clasp(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
    if (ep->dir != ep->prev->dir) return (EdgePtr)NULL;
    if (!perp[ep->dir][ep->next->dir]) return (EdgePtr)NULL;
    if (!anti[ep->prev->prev->dir][ep->next->dir]) return (EdgePtr)NULL;
    ivector loc;
    add_ivector(loc, ep->start, ep->next->increment);
    if (clkp->lattice[lat(loc)] == EMPTY) return (EdgePtr)NULL;
    ivector incr;
    copy_ivector(incr, increment_NEWSUD[kross[ep->dir][ep->next->dir]]);

    ivector loc1;

    add_ivector(loc1, ep->start, incr);
    if (clkp->lattice[lat(loc)] == EMPTY) return (EdgePtr)NULL;
    sub_ivector(loc1, ep->start, incr);
    if (clkp->lattice[lat(loc)] == EMPTY) return (EdgePtr)NULL;

    add_ivector(loc1, loc, incr);
    if (clkp->lattice[lat(loc)] == EMPTY) return (EdgePtr)NULL;
    sub_ivector(loc1, loc, incr);
    if (clkp->lattice[lat(loc)] == EMPTY) return (EdgePtr)NULL;

    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep1 = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (equal_ivector(loc, ep1->start))
            {
                if (ep1->dir != ep1->prev->dir) return (EdgePtr)NULL;
                if (!perp[ep1->dir][ep1->next->dir]) return (EdgePtr)NULL;
                if (!anti[ep1->prev->prev->dir][ep1->next->dir]) return (EdgePtr)NULL;
                ivector oloc;
                add_ivector(oloc, ep1->start, ep1->next->increment);
                if (!equal_ivector(oloc, ep->start)) return (EdgePtr)NULL;
                return ep1;
            }
            ep1 = ep1->next;
        }
        comp = comp->next;
    }
    return (EdgePtr)NULL; // should never get here
}

bool clk_check_path(char* path, CubicLatticeKnotPtr clkp, ivector loc, int value)
{
    ivector loco;
    int len = strlen(path);
    copy_ivector(loco, loc);
    for (int i = 0; i < len; i++)
    {
        add_ivector(loco, loco, increment_NEWSUD[(int)(path[i] - '0')]);
        if (clkp->lattice[lat(loco)] != value) return false;
    }
    return true;
}

EdgePtr is_parsite(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
    if (ep->dir != ep->prev->dir) return (EdgePtr)NULL;

    int xchange_dir;
    ivector loc;
    sub_ivector(loc, ep->start, clkp->loffset);

    switch (ep->dir)
    {
    case MOVE_NORTH:
    case MOVE_SOUTH:
        if (loc[0] > 0)
            xchange_dir = MOVE_EAST;
        else
            xchange_dir = MOVE_WEST;
        break;

    case MOVE_EAST:
    case MOVE_WEST:
        if (loc[1] > 0)
            xchange_dir = MOVE_NORTH;
        else
            xchange_dir = MOVE_SOUTH;
        break;

    default:
        return (EdgePtr)NULL;
    }


    char path[12];
    sprintf(path, "%d%d%d", xchange_dir, ep->dir, ep->dir);
    if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr)NULL;

    sprintf(path, "%d%d%d%d",
            MOVE_UP, ep->dir, xchange_dir, xchange_dir);
    if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr)NULL;

    ivector loc2;
    add_ivector(loc2, ep->start, increment_NEWSUD[xchange_dir]);
    add_ivector(loc2, loc2, increment_NEWSUD[xchange_dir]);
    sprintf(path, "%d%d%d%d",
            MOVE_DOWN, opposite[xchange_dir], opposite[xchange_dir], ep->dir);
    if (!clk_check_path(path, clkp, loc2, EMPTY)) return (EdgePtr)NULL;

    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep1 = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (equal_ivector(loc2, ep1->start))
            {
                if (ep1->dir != ep1->prev->dir || ep1->dir != ep->dir) return (EdgePtr)NULL;
                //printf ("dir found is %s\n", clk_dir_name [xchange_dir]); fflush (stdout);
                return ep1;
            }
            ep1 = ep1->next;
        }
        comp = comp->next;
    }
    return (EdgePtr)NULL; // should never get here
}

bool has_tight_clasp(CubicLatticeKnotPtr clkp)
{
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        EdgePtr epc;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (epc = is_tight_clasp(clkp, ep))
                return true;
            ep = ep->next;
        }
        comp = comp->next;
    }
    return false;
}

bool has_parsite(CubicLatticeKnotPtr clkp)
{
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        EdgePtr epc;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (epc = is_parsite(clkp, ep))
                return true;
            ep = ep->next;
        }
        comp = comp->next;
    }
    return false;
}

void collapse(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
    extern void delete_Edge(CubicLatticeKnotPtr knot, EdgePtr ep);
    //  printf ("collapsing edge #%d\n", ep->ID);
    clkp->lattice[lat(ep->prev->start)] = EMPTY;
    clkp->lattice[lat(ep->next->start)] = EMPTY;

    copy_ivector(ep->prev->start, ep->prev->prev->start);
    add_ivector(ep->start, ep->prev->start, ep->prev->increment);

    ComponentCLKPtr comp = (ComponentCLKPtr)ep->comp;
    delete_Edge(clkp, ep->prev->prev);
    delete_Edge(clkp, ep->next);
}

bool bfacf_strand_pass(CubicLatticeKnotPtr clkp)
{
    extern int rand_integer(int, int);
    if (!clkp) return false;
    if (clkp->nedges_total < 24) return false;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        EdgePtr epc;

        int shift = rand_integer(0, comp->nedges);
        for (int i = 0; i < shift; i++)
            ep = ep->next;

        for (int i = 0; i < comp->nedges; i++)
        {
            if (epc = is_tight_clasp(clkp, ep))
            {
                collapse(clkp, ep);
                collapse(clkp, epc);
                return true;
            }
            ep = ep->next;
        }
        comp = comp->next;
    }
    return false;
}

void join_edges(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
    extern bool clk_check_increment(ivector incr);
    extern int clk_direction(ivector incr);
    // assumes that start locations are valid
    ep1->next = ep2;
    ep2->prev = ep1;
    sub_ivector(ep1->increment, ep2->start, ep1->start);
    // bbbb
    if (!clk_check_increment(ep1->increment))
    {
        fprintf(stderr, "fails increment check!\n");
        exit(1);
    }
    ep1->dir = clk_direction(ep1->increment);
    clkp->lattice[lat(ep1->start)] = OCCUPIED;
}

void panic_exit_parsite(char* c)
{
    fprintf(stderr, "\n\n *** panick exit!  %s\n", c);
    exit(343);
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, EdgePtr ep)
{
    extern int rand_integer(int, int);
    int xchange_dir;
    EdgePtr epc;
    if (!(epc = is_parsite(clkp, ep))) return false;
    if (ep->dir != ep->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 1");
    if (ep->dir != epc->dir) panic_exit_parsite("bfacf_parsite_pass(): 2");
    if (epc->dir != epc->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 3");
    if (ep->dir == MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): 4");
    if (ep->dir == MOVE_DOWN) panic_exit_parsite("bfacf_parsite_pass(): 5");

    extern void edge_info(char* b, EdgePtr ep);
    edge_info("ep", ep);

    ivector diff;
    bool found = false;
    sub_ivector(diff, epc->start, ep->start);
    for (xchange_dir = 0; xchange_dir < 6; xchange_dir++)
    {
        ivector incr;
        mult_ivector(incr, increment_NEWSUD[xchange_dir], 2);
        if (equal_ivector(incr, diff))
        {
            found = true;
            break;
        }
    }
    // bbbb
    if (!found)
    {
        fprintf(stderr, "\n\n *** bfacf_parsite_pass(): xchange_dir not found!\n");
        fflush(stderr);
        exit(934);
    }

    EdgePtr epn = ep->next;
    EdgePtr epp = ep->prev;
    EdgePtr epcn = epc->next;
    EdgePtr epcp = epc->prev;

    // assign edges to nep[] array as follows:
    //          epp  ep   new edges  epn  epcp  epc  new edges   epcn
    // nep[]:     0   1       2 - 9   10    11   12    13 - 14      15

    // need 10 new edges
    // choose the first 10 unused edges in edge pool
    EdgePtr nep[16];
    int j = 2;
    for (int i = 0; i < 10; i++)
    {
        nep[j] = clkp->edgepool[clkp->nedges_total + i];
        nep[j]->comp = (void*)comp;
        j++;
        if (j == 10) j = 13;
    }
    nep[0] = epp;
    nep[1] = ep;
    nep[10] = epn;
    nep[11] = epcp;
    nep[12] = epc;
    nep[15] = epcn;

    int ep_orig_dir = ep->dir;
    // update positions
    add_ivector(nep[1]->start, nep[0]->start, increment_NEWSUD[MOVE_UP]);
    add_ivector(nep[2]->start, nep[1]->start, increment_NEWSUD[ep->dir]);
    add_ivector(nep[3]->start, nep[2]->start, increment_NEWSUD[xchange_dir]);
    add_ivector(nep[4]->start, nep[3]->start, increment_NEWSUD[xchange_dir]);
    add_ivector(nep[5]->start, nep[4]->start, increment_NEWSUD[MOVE_DOWN]);
    add_ivector(nep[6]->start, nep[5]->start, increment_NEWSUD[MOVE_DOWN]);
    add_ivector(nep[7]->start, nep[6]->start, increment_NEWSUD[opposite[xchange_dir]]);
    add_ivector(nep[8]->start, nep[7]->start, increment_NEWSUD[opposite[xchange_dir]]);
    add_ivector(nep[9]->start, nep[8]->start, increment_NEWSUD[MOVE_UP]);

    add_ivector(nep[12]->start, nep[11]->start, increment_NEWSUD[opposite[xchange_dir]]);
    add_ivector(nep[13]->start, nep[12]->start, increment_NEWSUD[ep->dir]);
    add_ivector(nep[14]->start, nep[13]->start, increment_NEWSUD[ep->dir]);

    // connect edges, adjust adjacency pointers and update lattice
    for (int i = 0; i < 15; i++)
    {
        if (i == 10) continue;
        join_edges(clkp, nep[i], nep[i + 1]);
    }

    // increase by 10 number of edges in comp and in knot
    comp->nedges += 10;
    clkp->nedges_total += 10;

    return true;
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp)
{
    extern int rand_integer(int, int);
    if (!clkp) return false;
    if (clkp->nedges_total < 24) return false;
    if (clkp->ncomps > 1) return false;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;

        int shift = rand_integer(0, comp->nedges);
        //printf ("shifting by %d\n", shift);
        for (int i = 0; i < shift; i++)
            ep = ep->next;

        for (int nn = 0; nn < comp->nedges; nn++)
        {
            if (bfacf_parsite_pass(clkp, comp, ep)) return true;
            ep = ep->next;
        }
        comp = comp->next;
    }
    return false;
}

int numb_tight_clasps(CubicLatticeKnotPtr clkp)
{
    return numb_tight_clasps(clkp, NULL);
}

int numb_tight_clasps(CubicLatticeKnotPtr clkp, void mark(ivector, ivector))
{
    int kount = 0;
    if (!clkp) return 0;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        EdgePtr epc;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (epc = is_tight_clasp(clkp, ep))
            {
                if (ep->ID > epc->ID)
                {
                    ++kount;
                    if (mark)
                        mark(ep->start, epc->start);
                }
            }
            ep = ep->next;
        }
        comp = comp->next;
    }
    return kount;
}

int numb_parsites(CubicLatticeKnotPtr clkp)
{
    return numb_parsites(clkp, NULL);
}

int numb_parsites(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr))
{
    int kount = 0;
    if (!clkp) return 0;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        EdgePtr epc;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (epc = is_parsite(clkp, ep))
            {
                if (ep->ID > epc->ID)
                {
                    ++kount;
                    if (mark)
                        mark(ep, epc);
                }
            }
            ep = ep->next;
        }
        comp = comp->next;
    }
    return kount;
}
