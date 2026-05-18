#include "bfacfMove.h"
#include "bfacf.h"
#include "clkConstants.h"

extern int rand_integer(int, int);
extern double rand_uniform();

// Lookup tables for BFACF moves

// given a direction DIR, opposite direction is opposite [DIR]
int opposite[6] = {MOVE_SOUTH, MOVE_WEST, MOVE_EAST, MOVE_NORTH, MOVE_DOWN, MOVE_UP};

ivector increment_NEWSUD[6] = {
    // increments associated with various directions
    {0, 1, 0}, // North
    {1, 0, 0}, // East
    {-1, 0, 0}, // West
    {0, -1, 0}, // South
    {0, 0, 1}, // Up
    {0, 0, -1}
}; // Down

int turn[6][4] = {
    // possible directions for +2 moves,
    // given direction of edge that is being moved
    {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // North
    {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // East
    {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // West
    {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // South
    {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}, // Up
    {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}
}; // Down

#define BAD  -1

int reflect[6][6] = {
    // first index is reflection direction, second index is direction to be reflected
    {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // North
    {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // East
    {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // West
    {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // South
    {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}, // Up
    {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}
}; // Down

int kross[6][6] = {
    {BAD, MOVE_DOWN, MOVE_UP, BAD, MOVE_EAST, MOVE_WEST}, // North
    {MOVE_UP, BAD, BAD, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
    {MOVE_DOWN, BAD, BAD, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
    {BAD, MOVE_UP, MOVE_DOWN, BAD, MOVE_WEST, MOVE_EAST}, // South
    {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, BAD, BAD}, // Up
    {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, BAD, BAD}
}; // Down

int rotate90[6][6] = {
    // first index is direction of axis of rotation, second index is direction to be rotated
    {MOVE_NORTH, MOVE_DOWN, MOVE_UP, MOVE_SOUTH, MOVE_EAST, MOVE_WEST}, // North
    {MOVE_UP, MOVE_EAST, MOVE_WEST, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
    {MOVE_DOWN, MOVE_EAST, MOVE_WEST, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
    {MOVE_NORTH, MOVE_UP, MOVE_DOWN, MOVE_SOUTH, MOVE_WEST, MOVE_EAST}, // South
    {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, MOVE_UP, MOVE_DOWN}, // Up
    {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, MOVE_UP, MOVE_DOWN}
}; // Down

bool perp[6][6] = {
    // given direction DIRA and DIRB, perp [DIRA][DIRB] is true
    // iff DIRA perpendicular to DIRB
    {false, true, true, false, true, true}, // North
    {true, false, false, true, true, true}, // East
    {true, false, false, true, true, true}, // West
    {false, true, true, false, true, true}, // South
    {true, true, true, true, false, false}, // Up
    {true, true, true, true, false, false}
}; // Down

bool anti[6][6] = {
    // given direction DIRA and DIRB, anti [DIRA][DIRB] is true
    // iff DIRA antiparallel to DIRB
    {false, false, false, true, false, false}, // North
    {false, false, true, false, false, false}, // East
    {false, true, false, false, false, false}, // West
    {true, false, false, false, false, false}, // South
    {false, false, false, false, false, true}, // Up
    {false, false, false, false, true, false}
}; // Down

// names of the directions
char clk_dir_name[6][6] = {"North", "East", "West", "South", "Up", "Down"};

void delete_Edge(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
    // adjust adjacent edges' adjacency pointers
    ep->prev->next = ep->next;
    ep->next->prev = ep->prev;

    // swap this edge with the last edge in the edgepool
    EdgePtr tmp = knot->edgepool[knot->nedges_total - 1];
    tmp->locpool = ep->locpool;
    knot->edgepool[ep->locpool] = tmp;
    knot->edgepool[knot->nedges_total - 1] = ep;
    ep->locpool = knot->nedges_total - 1;

    // reduce by 1 the number of edges in knot and in comp
    --knot->nedges_total;
    --comp->nedges;

    // update info contained in `comp', in the event we've
    // deleted first or last edge in that component

    if (comp->first_edge == ep)
        comp->first_edge = ep->next;
    else if (comp->last_edge == ep)
        comp->last_edge = ep->prev;
}

bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep, ivector increment, int dir)
{
    // increment is the direction of the move
    // dir is the direction associated with this increment

    if (knot->nedges_total > knot->poolsize - 3)
    {
        return false;
    }
    if (comp->nedges > comp->maxedges) return false;

    // location that edge might want to move to must be empty

    ivector test_locationA, test_locationB;
    add_ivector(test_locationA, ep->start, increment);

    if (clk_check_for_edge_hits(knot, dir, test_locationA))
    {
        if (knot->auto_recentre)
        {
            if (recentre_knot_in_lattice(knot)) // returns true if knot can't be recentered
                return false;
            add_ivector(test_locationA, ep->start, increment); // need to recalculate this
        }
        else
            return false;
    }

    if (knot->lattice[lat(test_locationA)] == OCCUPIED) return false;
    add_ivector(test_locationB, ep->next->start, increment);
    if (knot->lattice[lat(test_locationB)] == OCCUPIED) return false;

    // choose the first two unused edges in edge pool
    EdgePtr ep_prev = knot->edgepool[knot->nedges_total]; // note: not nedges - 1 as in above
    EdgePtr ep_next = knot->edgepool[knot->nedges_total + 1]; // note: not nedges - 1 as in above
    ep_prev->comp = ep_next->comp = (void*)comp;

    // adjust adjacency pointers
    ep_prev->prev = ep->prev;
    ep_prev->next = ep;
    ep_next->prev = ep;
    ep_next->next = ep->next;
    ep->prev->next = ep_prev;
    ep->next->prev = ep_next;
    ep->prev = ep_prev;
    ep->next = ep_next;

    copy_ivector(ep_prev->increment, increment);
    negate_ivector(ep_next->increment, increment);
    ep_prev->dir = dir;
    ep_next->dir = opposite[dir];


    // increase by 2 number of edges in comp and in knot
    comp->nedges += 2;
    knot->nedges_total += 2;

    // update positions

    EdgePtr ep2 = ep_prev;
    for (int i = 0; i < 3; i++)
    {
        add_ivector(ep2->start, ep2->prev->start, ep2->prev->increment);
        ep2 = ep2->next;
    }

    // update lattice

    knot->lattice[lat(ep->start)] = OCCUPIED;
    knot->lattice[lat(ep->next->start)] = OCCUPIED;


    knot->success_plus2++;
#ifdef TESTING
    knot->p2dir[dir]++;
#endif
    return true;
}

bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
    // this function assumes that a move in any of the four directions
    // is possible locally (ignoring global self-intersections)

    // choose a direction at random
    int dir = turn[ep->dir][rand_integer(0, 4)];
    return perform_plus2_move(knot, comp, ep, increment_NEWSUD[dir], dir);
}

bool perform_plus2_move_alt(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
    int dir;
    // this function is similar to perform_plus2_move() above, but
    // assumes that a move in at most three of the four directions
    // is possible locally (ignoring global self-intersections)

    // in order for the move to be possible locally, direction of
    // the move cannot be in the same direction as the next edge,
    // and not in the opposite direction as the previous edge

    do
    {
        dir = turn[ep->dir][rand_integer(0, 4)];
    }
    while (dir == ep->next->dir || dir == opposite[ep->prev->dir]);

    return perform_plus2_move(knot, comp, ep, increment_NEWSUD[dir], dir);
}

bool perform_0_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
    // type 0 moves are only possible if this edge and the next edge
    // are perpendicular

    if (ep->frozen || ep->next->frozen) return false;

    // also the new lattice location needs to be empty
    ivector test_location;
    add_ivector(test_location, ep->start, ep->next->increment);
    if (knot->lattice[lat(test_location)] == OCCUPIED) return false;

    ivector tmp;
    // simply swap increments between this edge and next
    copy_ivector(tmp, ep->increment);
    copy_ivector(ep->increment, ep->next->increment);
    copy_ivector(ep->next->increment, tmp);

    // clear lattice at old vertex location
    knot->lattice[lat(ep->next->start)] = EMPTY;

    int itmp;
    itmp = ep->dir;
    ep->dir = ep->next->dir;
    ep->next->dir = itmp;

    // update position and lattice
    add_ivector(ep->next->start, ep->start, ep->increment);
    knot->lattice[lat(ep->next->start)] = OCCUPIED;

    knot->success_0++;

    return true;
}

bool perform_minus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
    if (comp->nedges < 5) return false;
    if (comp->nedges < comp->minedges) return false;
    if (ep->next->frozen || ep->prev->frozen) return false;

    // clear lattice
    knot->lattice[lat(ep->start)] = EMPTY;
    knot->lattice[lat(ep->next->start)] = EMPTY;

    // for this case, no increments need to be changed or specified
    delete_Edge(knot, comp, ep->prev);
    delete_Edge(knot, comp, ep->next);

    // update position
    add_ivector(ep->start, ep->prev->start, ep->prev->increment);

    knot->success_minus2++;
    return true;
}

bool perform_move(CubicLatticeKnotPtr knot)
{
    if (knot->nfrozen == knot->nedges_total) return false;

    // choose an edge uniformly at random from the edge pool

    EdgePtr ep = knot->edgepool[rand_integer(knot->nfrozen, knot->nedges_total)];
    ComponentCLKPtr comp = (ComponentCLKPtr)ep->comp;

#ifdef CHECK_CONFIG_not_used
    ivector eploc;
    copy_ivector(eploc, ep->start);
#endif

    // The types of moves possible depend on the chosen edge, ep, and
    // the two adjacent edges, ep->prev and ep->next.

    double p = rand_uniform(); // uniform number between 0 and 1
    bool value = false; // assume move fails

    // NOTE: reordering the cases below may increase performance slightly
    // but Case 1 needs to be tested before Case 3.

    // Case 1: chosen edge is parallel to both adjacent edges,
    // all four moves are +2 moves.

    // directions of adjacent edges same as chosen edge

    if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
    {
        if (p < comp->p_4p2)
            value = perform_plus2_move(knot, comp, ep);
    }

    // Case 2: chosen edge is perpendicular to both adjacent edges
    // and those adjacent edges are anti-parallel to each other.
    // One move is a -2 move, other three are +2 moves.

    else if (perp[ep->dir][ep->next->dir] && anti[ep->prev->dir][ep->next->dir])
    {
        if (p < comp->p_minus2)
            value = perform_minus2_move(knot, comp, ep);
        else if (p < comp->p_m23p2)
            value = perform_plus2_move_alt(knot, comp, ep);
    }

    // Case 3: chosen edge is perpendicular to one adjacent edge
    // and parallel to the other adjacent edge.
    // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
    // because the case of both being parallel has been ruled out above.
    // One move is a 0 move, other three are +2 moves.

    // NOTE: test this before Case 1 above

    else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
    {
        if (p < comp->p_0)
        {
            if (ep->dir == ep->prev->dir)
                value = perform_0_move(knot, comp, ep);
            else
                value = perform_0_move(knot, comp, ep->prev);
        }
        else if (p < comp->p_03p2)
            value = perform_plus2_move_alt(knot, comp, ep);
    }

    // Case 4: Only possibility remaining: chosen edge is perpendicular
    // to both adjacent edges, and the adjacent edges are not anti-parallel.
    // Two moves are 0 moves and two moves are +2 moves.

    else
    {
        if (p < comp->p_2p0)
        {
            if (p < comp->p_0)
                value = perform_0_move(knot, comp, ep->prev);
            else
                value = perform_0_move(knot, comp, ep);
        }
        else if (p < comp->p_2p02p2)
            value = perform_plus2_move_alt(knot, comp, ep);
    }

#ifdef CHECK_CONFIG_NOT_USED
    extern void bfacf_check_config(ivector);
    bfacf_check_config(eploc);
#endif

#ifdef TESTING2    // Warning: this changes the algorithm from roughly constant time
    //          to roughly linear time (in number of edges)
    if (value)
        clk_check_increments(knot);

#endif

    return value;
}
