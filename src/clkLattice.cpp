#include "clkLattice.h"
#include "clkConstants.h"
#include "clkTables.h"


void init_lattice(CubicLatticeKnotPtr clkp, ivector mid)
{
    // offset the knot so that the point `mid' lies at at
    // (LATTICE_SIZE/2, LATTICE_SIZE/2, LATTICE_SIZE/2)
    set_ivector(clkp->loffset, LATTICE_SIZE / 2, LATTICE_SIZE / 2, LATTICE_SIZE / 2);
    sub_ivector(clkp->loffset, clkp->loffset, mid);
    set_ivector(clkp->max_range, LATTICE_SIZE - 2, LATTICE_SIZE - 2, LATTICE_SIZE - 2);
    ComponentCLKPtr comp = clkp->fcomp;

    while (comp)
    {
        // now go thru all the edges, adding the loffset and marking as OCCUPIED
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            add_ivector(ep->start, ep->start, clkp->loffset); // we will subtract this offset when saving the knot
            clkp->lattice[lat(ep->start)] = OCCUPIED;
            ep = ep->next;
        }
        comp = comp->next;
    }
}

bool recentre_knot_in_lattice(CubicLatticeKnotPtr knot)
{
    // returns true if knot can't be recentered
    if (!knot) return true;

    ComponentCLKPtr comp;
    ivector min, max, mid;
    ivector diff;
    EdgePtr ep;

    extern void clk_get_extent(ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot);
    clk_get_extent(min, max, knot);
    if (max[0] - min[0] > knot->max_range[0]) return true;
    if (max[1] - min[1] > knot->max_range[1]) return true;
    if (max[2] - min[2] > knot->max_range[2]) return true;

    // first clear the old path
    clear_lattice(knot);
    ep = knot->fcomp->first_edge;

#ifdef RGS_PARANOIA
    sub_ivector(diff, ep->start, knot->loffset);
    printf("before: %d %d %d   %d %d %d   %d %d %d\n",
           ep->start[0], ep->start[1], ep->start[2],
           knot->loffset[0], knot->loffset[1], knot->loffset[2],
           diff[0], diff[1], diff[2]);
    fflush(stdout);
#endif

    comp = knot->fcomp;
    while (comp)
    {
        ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            // recover the start locations for the edges in actual 3D space
            sub_ivector(ep->start, ep->start, knot->loffset);
            ep = ep->next;
        }
        comp = comp->next;
    }

    extern void midpoint(ivector res, ivector a, ivector b);
    midpoint(mid, max, min);

    init_lattice(knot, mid);

#ifdef RGS_PARANOIA
    ep = knot->fcomp->first_edge;
    sub_ivector(diff, ep->start, knot->loffset);
    printf("after: %d %d %d   %d %d %d   %d %d %d\n\n",
           ep->start[0], ep->start[1], ep->start[2],
           knot->loffset[0], knot->loffset[1], knot->loffset[2],
           diff[0], diff[1], diff[2]);
    fflush(stdout);
#endif

    return false;
}

void set_lattice(CubicLatticeKnotPtr clkp, char value, ivector min, ivector max)
{
    if (!clkp) return;
    //clkp->auto_recentre = false;
    ivector offset_min, offset_max;
    add_ivector(offset_min, min, clkp->loffset);
    add_ivector(offset_max, max, clkp->loffset);
    for (int i = 0; i < 3; i++)
    {
        CLAMP(offset_min [i], 0, LATTICE_SIZE - 1);
        CLAMP(offset_max [i], 0, LATTICE_SIZE - 1);
    }
    ivector a;

    for (a[0] = offset_min[0]; a[0] <= offset_max[0]; a[0]++)
        for (a[1] = offset_min[1]; a[1] <= offset_max[1]; a[1]++)
            for (a[2] = offset_min[2]; a[2] <= offset_max[2]; a[2]++)
                clkp->lattice[lat(a)] = value;
}

void fill_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max)
{
    set_lattice(clkp, (char)OCCUPIED, min, max);
}

void invert_lattice(CubicLatticeKnotPtr clkp)
{
    if (!clkp) return;
    for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
    {
        if (clkp->lattice[i] == OCCUPIED)
            clkp->lattice[i] = EMPTY;
        else
            clkp->lattice[i] = OCCUPIED;
    }
}

void dilate_lattice(CubicLatticeKnotPtr clkp)
{
    // expand occupied region to occupy neighboring cells
    if (!clkp) return;
    ivector pos;

    for (pos[2] = 1; pos[2] < LATTICE_SIZE - 1; pos[2]++)
    {
        // unnecessesary brace to make VC++ indenter happy
        for (pos[1] = 1; pos[1] < LATTICE_SIZE - 1; pos[1]++)
        {
            for (pos[0] = 1; pos[0] < LATTICE_SIZE - 1; pos[0]++)
            {
                if (clkp->lattice[lat(pos)] == OCCUPIED)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        ivector pos2;
                        add_ivector(pos2, pos, increment_NEWSUD[i]);
                        if (clkp->lattice[lat(pos2)] == EMPTY)
                        {
                            clkp->lattice[lat(pos2)] = NEXTOCCUPIED;
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
    {
        if (clkp->lattice[i] == NEXTOCCUPIED)
            clkp->lattice[i] = OCCUPIED;
    }
}

int bfacf_info_filled(CubicLatticeKnotPtr knot)
{
    int filled = 0;
    for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
        if (knot->lattice[i] == OCCUPIED)
            ++filled;
    return filled;
}

void set_lattice(CubicLatticeKnotPtr clkp, char value)
{
    if (!clkp) return;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        ivector a;
        for (int i = 0; i < comp->nedges; i++)
        {
            clkp->lattice[lat(ep->start)] = value;
            ep = ep->next;
        }
        if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
        {
            add_ivector(a, ep->start, ep->increment);
            clkp->lattice[lat(a)] = value;
        }
        comp = comp->next;
    }
}

void fill_lattice(CubicLatticeKnotPtr clkp)
{
    set_lattice(clkp, OCCUPIED);
}

void clear_lattice(CubicLatticeKnotPtr clkp)
{
    set_lattice(clkp, EMPTY);
}

void clear_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max)
{
    set_lattice(clkp, (char)EMPTY, min, max);
}
