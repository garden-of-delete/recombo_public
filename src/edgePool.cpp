#include "edgePool.h"
#include "clkConstants.h"
#include "clkTables.h"
#include "clkGeometry.h"

int clk_minedges = 0;
int clk_maxedges = HUGE_NUMBER;

void swapEdgesPool(CubicLatticeKnotPtr clkp, int loc1, int loc2)
{
    if (loc1 == loc2) return;
    EdgePtr ept1 = clkp->edgepool[loc1];
    EdgePtr ept2 = clkp->edgepool[loc2];
    clkp->edgepool[loc1] = ept2;
    clkp->edgepool[loc2] = ept1;
    ept1->locpool = loc2;
    ept2->locpool = loc1;
}

void freezeEdge(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
    if (ep->locpool < clkp->nfrozen) // edge is already frozen
        return;
    swapEdgesPool(clkp, ep->locpool, clkp->nfrozen);
    clkp->nfrozen++;
}

bool freezeEdge(CubicLatticeKnotPtr clkp, ivector start)
{
    // freeze edge in `clkp' that starts at `start'
    if (!clkp) return false;
    EdgePtr ep = clk_get_edge(clkp, start);
    if (!ep) return false;
    freezeEdge(clkp, ep);
    ep->frozen = true;
    return true;
}

void freezeEdge(EdgePtr ep)
{
    ComponentCLKPtr comp = (ComponentCLKPtr)ep->comp;
    CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr)comp->clkp;
    freezeEdge(clkp, ep);
}

void thawEdge(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
    if (ep->locpool >= clkp->nfrozen) // edge is already thawed
        return;
    swapEdgesPool(clkp, ep->locpool, clkp->nfrozen - 1);
    clkp->nfrozen--;
}

void thawEdge(EdgePtr ep)
{
    ComponentCLKPtr comp = (ComponentCLKPtr)ep->comp;
    CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr)comp->clkp;
    thawEdge(clkp, ep);
}

bool freeze_component(CubicLatticeKnotPtr clkp, int ID, bool freeze)
{
    if (!clkp) return false;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        if (ID == comp->ID)
        {
            EdgePtr ep = comp->first_edge;
            for (int i = 0; i < comp->nedges; i++)
            {
                if (freeze)
                    freezeEdge(clkp, ep);
                else
                    thawEdge(clkp, ep);
                ep = ep->next;
            }
            return true;
        }
        comp = comp->next;
    }
    return false;
}

void clk_set_edge_limit(int min, int max)
{
    clk_minedges = min;
    clk_maxedges = max;
}

void clk_set_nedge_limit(CubicLatticeKnotPtr clkp, int min, int max, int ID)
{
    if (!clkp) return;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        if (ID < 0 || ID == comp->ID)
        {
            comp->minedges = min;
            comp->maxedges = max;
        }
        comp = comp->next;
    }
}
