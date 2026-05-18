#include "clkGeometry.h"
#include "clkConstants.h"

EdgePtr clk_get_edge(CubicLatticeKnotPtr clkp, ivector start)
{
    ComponentCLKPtr comp = clkp->fcomp;
    add_ivector(start, start, clkp->loffset);
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            if (equal_ivector(ep->start, start))
            {
                return ep;
            }
            ep = ep->next;
        }
        comp = comp->next;
    }
    return (EdgePtr)NULL;
}

int clk_count_edges(EdgePtr ep)
{
    int kount = 1;
    EdgePtr start = ep;
    ep = ep->next;
    while (ep != start)
    {
        kount++;
        ep = ep->next;
    }
    return kount;
}

int clk_arc_length_distance(EdgePtr ep1, EdgePtr ep2)
{
    ComponentCLKPtr comp1 = (ComponentCLKPtr)ep1->comp;
    ComponentCLKPtr comp2 = (ComponentCLKPtr)ep2->comp;
    if (comp1 != comp2)
        return ARC_LENGTH_DISTANCE_INFINITE;

    EdgePtr epn = ep1->next;
    EdgePtr epp = ep1->prev;
    int distance = 1;
    while (epn != ep2 && epp != ep2)
    {
        epn = epn->next;
        epp = epp->prev;
        ++distance;
    }

    return distance;
}

int clk_arc_length_distance(CubicLatticeKnotPtr clkp, ivector s1, ivector s2)
{
    EdgePtr ep1 = clk_get_edge(clkp, s1);
    if (!ep1) return ARC_LENGTH_DISTANCE_NO_EDGE1;
    EdgePtr ep2 = clk_get_edge(clkp, s2);
    if (!ep2) return ARC_LENGTH_DISTANCE_NO_EDGE2;
    return clk_arc_length_distance(ep1, ep2);
}

void clk_get_extent(ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot)
{
    if (!knot) return;
    copy_ivector(minbb, knot->fcomp->first_edge->start);
    copy_ivector(maxbb, minbb);
    ComponentCLKPtr comp = knot->fcomp;

    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        for (int i = 0; i < comp->nedges; i++)
        {
            min_ivector(minbb, minbb, ep->start);
            max_ivector(maxbb, maxbb, ep->start);
            ep = ep->next;
        }
        comp = comp->next;
    }
    sub_ivector(minbb, minbb, knot->loffset);
    sub_ivector(maxbb, maxbb, knot->loffset);
}

void clk_bounding_box(int& X, int& Y, int& Z, CubicLatticeKnotPtr knot)
{
    ivector minbb, maxbb;
    clk_get_extent(minbb, maxbb, knot);
    X = maxbb[0] - minbb[0];
    Y = maxbb[1] - minbb[1];
    Z = maxbb[2] - minbb[2];
}

int clk_direction(ivector incr)
{
    // return the direction associated with an increment
#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

    switch ((1 << (incr[0] + 1)) | (1 << (incr[1] + 1 + 2)) | (1 << (incr[2] + 1 + 4)))
    {
    case NORTH:
        return MOVE_NORTH;
        break;
    case EAST:
        return MOVE_EAST;
        break;
    case WEST:
        return MOVE_WEST;
        break;
    case SOUTH:
        return MOVE_SOUTH;
        break;
    case UP:
        return MOVE_UP;
        break;
    case DOWN:
        return MOVE_DOWN;
        break;
    }
    return MOVE_INVALID;
}

int clk_direction(ivector end, ivector start)
{
    ivector incr;
    sub_ivector(incr, end, start);
    return clk_direction(incr);
}
