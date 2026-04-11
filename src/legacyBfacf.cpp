#include "legacyBfacf.h"

#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

// from bfacf.cpp

void delete_Edge(CubicLatticeKnotPtr knot, EdgePtr ep);

// directions are specified by an integer from 0 to 5 inclusive
// MOVE_NORTH is 0, MOVE_EAST is 1, etc.  (see bfacf.h for definitions)

// given a direction DIR, opposite direction is opposite [DIR]
int opposite [6] = {MOVE_SOUTH, MOVE_WEST, MOVE_EAST, MOVE_NORTH, MOVE_DOWN, MOVE_UP};

ivector increment_NEWSUD [6] = {// increments associated with various directions
   { 0, 1, 0}, // North
   { 1, 0, 0}, // East
   {-1, 0, 0}, // West
   { 0, -1, 0}, // South
   { 0, 0, 1}, // Up
   { 0, 0, -1}
}; // Down

int turn [6][4] = {// possible directions for +2 moves, 
   // given direction of edge that is being moved 
   {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // North
   {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // East
   {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // West
   {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // South
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}, // Up
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}
}; // Down


#define BAD  -1

int reflect [6][6] = {// first index is reflection direction, second index is direction to be reflected
   {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // North
   {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // East
   {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // West
   {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // South
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}, // Up
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}
}; // Down

int kross [6][6] = {
   {BAD, MOVE_DOWN, MOVE_UP, BAD, MOVE_EAST, MOVE_WEST}, // North
   {MOVE_UP, BAD, BAD, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
   {MOVE_DOWN, BAD, BAD, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
   {BAD, MOVE_UP, MOVE_DOWN, BAD, MOVE_WEST, MOVE_EAST}, // South
   {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, BAD, BAD}, // Up
   {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, BAD, BAD}
}; // Down

int rotate90 [6][6] = {// first index is direction of axis of rotation, second index is direction to be rotated
   {MOVE_NORTH, MOVE_DOWN, MOVE_UP, MOVE_SOUTH, MOVE_EAST, MOVE_WEST}, // North
   {MOVE_UP, MOVE_EAST, MOVE_WEST, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
   {MOVE_DOWN, MOVE_EAST, MOVE_WEST, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
   {MOVE_NORTH, MOVE_UP, MOVE_DOWN, MOVE_SOUTH, MOVE_WEST, MOVE_EAST}, // South
   {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, MOVE_UP, MOVE_DOWN}, // Up
   {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, MOVE_UP, MOVE_DOWN}
}; // Down

bool perp [6][6] = {// given direction DIRA and DIRB, perp [DIRA][DIRB] is true 
   // iff DIRA perpendicular to DIRB
   {false, true, true, false, true, true}, // North
   {true, false, false, true, true, true}, // East
   {true, false, false, true, true, true}, // West
   {false, true, true, false, true, true}, // South
   {true, true, true, true, false, false}, // Up
   {true, true, true, true, false, false}
}; // Down

bool anti [6][6] = {// given direction DIRA and DIRB, anti [DIRA][DIRB] is true 
   // iff DIRA antiparallel to DIRB
   {false, false, false, true, false, false}, // North
   {false, false, true, false, false, false}, // East
   {false, true, false, false, false, false}, // West
   {true, false, false, false, false, false}, // South
   {false, false, false, false, false, true}, // Up
   {false, false, false, false, true, false}
}; // Down

// names of the directions
char clk_dir_name [6][6] = {"North", "East", "West", "South", "Up", "Down"};

bool add_edge_to_knot(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end)
{
   EdgePtr ep = (EdgePtr) calloc(1, sizeof (Edge)); // create a new edge

   sub_ivector(ep->increment, end, start); // increment is (end - start)
   ep->dir = clk_direction(ep->increment); // direction associated with increment
   copy_ivector(ep->start, start); // starting location of this edge
   ep->frozen = false; // edge is free to move
   ep->comp = (void *) comp;

   if (comp->first_edge)
   { // some edges have already been added
      comp->last_edge->next = ep;
      ep->prev = comp->last_edge;
      comp->last_edge = ep;
   }
   else
   { // this is the first edge
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
   clkp->edgepool = (EdgePtr *) calloc(clkp->poolsize, sizeof (EdgePtr));
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
         clkp->edgepool [index++] = ep;
         ep = ep->next;
      }
      comp = comp->next;
   }

   while (index < poolsize)
      clkp->edgepool [index++] = (EdgePtr) calloc(1, sizeof (Edge)); // new, unused edge

   for (index = 0; index < clkp->poolsize; index++)
   {
      clkp->edgepool [index]->locpool = index;
      clkp->edgepool [index]->ID = index;
   }

   clkp->lattice = (char *) calloc(LATTICE_TOTAL_SIZE, sizeof (char));
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

// the three moves below assume that that move is legal locally, but maybe not globally

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

   if (knot->lattice [lat(test_locationA)] == OCCUPIED) return false;
   add_ivector(test_locationB, ep->next->start, increment);
   if (knot->lattice [lat(test_locationB)] == OCCUPIED) return false;

   // choose the first two unused edges in edge pool
   EdgePtr ep_prev = knot->edgepool [knot->nedges_total]; // note: not nedges - 1 as in above
   EdgePtr ep_next = knot->edgepool [knot->nedges_total + 1]; // note: not nedges - 1 as in above
   ep_prev->comp = ep_next->comp = (void *) comp;

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
   ep_next->dir = opposite [dir];


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

   knot->lattice [lat(ep->start)] = OCCUPIED;
   knot->lattice [lat(ep->next->start)] = OCCUPIED;


   knot->success_plus2++;
#ifdef TESTING
   knot->p2dir [dir]++;
#endif
   return true;
}

bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // this function assumes that a move in any of the four directions
   // is possible locally (ignoring global self-intersections)

   // choose a direction at random
   int dir = turn [ep->dir][rand_integer(0, 4)];
   return perform_plus2_move(knot, comp, ep, increment_NEWSUD [dir], dir);
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
      dir = turn [ep->dir][rand_integer(0, 4)];
   }
   while (dir == ep->next->dir || dir == opposite [ep->prev->dir]);

   return perform_plus2_move(knot, comp, ep, increment_NEWSUD [dir], dir);
}

bool perform_0_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // type 0 moves are only possible if this edge and the next edge
   // are perpendicular

   if (ep->frozen || ep->next->frozen) return false; // NEW!!

   // also the new lattice location needs to be empty
   ivector test_location;
   add_ivector(test_location, ep->start, ep->next->increment);
   if (knot->lattice [lat(test_location)] == OCCUPIED) return false;

   ivector tmp;
   // simply swap increments between this edge and next
   copy_ivector(tmp, ep->increment);
   copy_ivector(ep->increment, ep->next->increment);
   copy_ivector(ep->next->increment, tmp);

   // clear lattice at old vertex location
   knot->lattice [lat(ep->next->start)] = EMPTY;

   int itmp;
   itmp = ep->dir;
   ep->dir = ep->next->dir;
   ep->next->dir = itmp;

   // update position and lattice
   add_ivector(ep->next->start, ep->start, ep->increment);
   knot->lattice [lat(ep->next->start)] = OCCUPIED;

   knot->success_0++;

   return true;
}

void delete_Edge(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // adjust adjacent edges' adjacency pointers
   ep->prev->next = ep->next;
   ep->next->prev = ep->prev;

   // swap this edge with the last edge in the edgepool
   EdgePtr tmp = knot->edgepool [knot->nedges_total - 1];
   tmp->locpool = ep->locpool;
   knot->edgepool [ep->locpool] = tmp;
   knot->edgepool [knot->nedges_total - 1] = ep;
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

bool perform_minus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   if (comp->nedges < 5) return false;
   if (comp->nedges < comp->minedges) return false;
   if (ep->next->frozen || ep->prev->frozen) return false;

   // clear lattice
   knot->lattice [lat(ep->start)] = EMPTY;
   knot->lattice [lat(ep->next->start)] = EMPTY;

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

   EdgePtr ep = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;

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

   else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
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
   return (EdgePtr) NULL;
}

void swapEdgesPool(CubicLatticeKnotPtr clkp, int loc1, int loc2)
{
   if (loc1 == loc2) return;
   EdgePtr ept1 = clkp->edgepool [loc1];
   EdgePtr ept2 = clkp->edgepool [loc2];
   clkp->edgepool [loc1] = ept2;
   clkp->edgepool [loc2] = ept1;
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
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr) comp->clkp;
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
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr) comp->clkp;
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

int clk_minedges = 0;
int clk_maxedges = HUGE_NUMBER;

void clk_set_edge_limit(int min, int max)
{ // for newly created components
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

// from clk_util.cpp

bool clk_recombo_limit = false;
int clk_min_arc_length_distance = 0;

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

int clk_validate(CubicLatticeKnotPtr clkp) {
    for (int i = 0; i < clkp->nedges_total; i++)
        if (clkp->edgepool[i]->locpool != i)
            return 1;

   int kount = 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr start;
      EdgePtr ep = start = comp->first_edge;
      if (comp->nedges < 4) return 7;
      ivector loc;
      copy_ivector(loc, ep->start);
      for (int i = 0; i < comp->nedges; i++)
      {
         if (ep->prev->next != ep)
            return 2;
         if (ep->next->prev != ep)
            return 3;
         if ((ComponentCLKPtr) ep->comp != comp)
            return 5;
         if (clkp->lattice [lat(ep->start)] != OCCUPIED)
            return 8;
         add_ivector(loc, loc, ep->increment);
         if (!equal_ivector(loc, ep->next->start))
            return 9;
         if (ep->dir != clk_direction(ep->increment))
            return 10;
         ep = ep->next;
      }
      if (ep != start) return 4;
      kount += comp->nedges;
      comp = comp->next;
   }
   if (kount != clkp->nedges_total) return 6;
   return 0;
}

void clk_validate(CubicLatticeKnotPtr clkp, char *s)
{
   int code = clk_validate(clkp);
   if (code)
   {
      fprintf(stderr, "fails at %s with code %d\n", s, code);
      exit(102);
   }
   //  else 
   //  fprintf (stderr, "ok at %s\n", s);
}

bool clk_check_increment(ivector incr)
{
   // a valid increment should be +-1 in one direction only
   if (abs(incr [0]) + abs(incr [1]) + abs(incr [2]) == 1)
      return true;
   else
      return false;
}

bool clk_check_lattice(CubicLatticeKnotPtr knot)
{
   if (!knot) return false;
   // WARNING: this function has the side effect of clearing the knot path!
   //          use only at the end of a run to make sure things are OK

   ComponentCLKPtr comp = knot->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      printf("checking lattice ... ");
      fflush(stdout);
      for (int e = 0; e < comp->nedges; e++)
      {
         if (knot->lattice [lat(ep->start)] != OCCUPIED)
         {
            printf("\n\nlattice location %d %d %d for edge with ID %d is not occupied!\n",
                    ep->start [0], ep->start [1], ep->start [2], ep->ID);
            fflush(stdout);
            exit(102);
         }
         //printf ("%d ", ep->ID); fflush (stdout);
         knot->lattice [lat(ep->start)] = EMPTY;

         ep = ep->next;
      }

      if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
      {
         ivector endpoint;
         add_ivector(endpoint, comp->last_edge->start, comp->last_edge->increment);
         knot->lattice [lat(endpoint)] = EMPTY;
      }

      comp = comp->next;
   }

   int ix, iy, iz;
   ivector loc;

   // we've cleared the knot path, so every location in lattice should now be EMPTY

   for (ix = 0; ix < LATTICE_SIZE; ix++)
   { // unnecessesary brace to make VC++ indenter happy
      for (iy = 0; iy < LATTICE_SIZE; iy++)
      {
         for (iz = 0; iz < LATTICE_SIZE; iz++)
         {
            set_ivector(loc, ix, iy, iz);
            if (knot->lattice [lat(loc)] != EMPTY)
            {
               printf("\n\nlattice location %d %d %d is not EMPTY as it should be!\n\n",
                       ix, iy, iz);
               fflush(stdout);
               exit(102);
            }
         }
      }
   }

   printf("lattice checks out!\n");
   return true;
}

int clk_arc_length_distance(EdgePtr ep1, EdgePtr ep2)
{
   ComponentCLKPtr comp1 = (ComponentCLKPtr) ep1->comp;
   ComponentCLKPtr comp2 = (ComponentCLKPtr) ep2->comp;
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

void clk_bounding_box(int &X, int &Y, int &Z, CubicLatticeKnotPtr knot)
{
   ivector minbb, maxbb;
   clk_get_extent(minbb, maxbb, knot);
   X = maxbb [0] - minbb [0];
   Y = maxbb [1] - minbb [1];
   Z = maxbb [2] - minbb [2];
}

void clk_check_increments(CubicLatticeKnotPtr knot)
{
   if (knot->nedges_total < 4)
   {
      fprintf(stderr,
              "\n\n*** clk_check_increments(): Ridiculous number of edges (%d) in knot!\n\n",
              knot->nedges_total);
      exit(102);
   }

   ComponentCLKPtr comp = knot->fcomp;
   FILE *fp = (FILE *) NULL;
   //fp = fopen ("out.txt", "w");

   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      ivector pos, start;

      copy_ivector(pos, comp->first_edge->start);
      copy_ivector(start, pos);

      int N = comp->nedges;
      if (comp->flags & COMPONENT_CLK_FLAG_OPEN) N--;

      for (int i = 0; i < N; i++)
      {
         if (!clk_check_increment(ep->increment))
         {
            fprintf(stderr, "\n\n clk_check_increments(): bad increment of (%d, %d, %d) seen!\n\n",
                    ep->increment [0], ep->increment [1], ep->increment [2]);
            exit(343);
         }
         if (clk_direction(ep->increment) != ep->dir)
         {
            fprintf(stderr, "\n\n clk_check_increments(): increment of (%d, %d, %d) is not direction %s!!\n\n",
                    ep->increment [0], ep->increment [1], ep->increment [2], clk_dir_name [ep->dir]);
            exit(93454);
         }
         add_ivector(pos, pos, ep->increment);
         if (!equal_ivector(pos, ep->next->start))
         {
            fprintf(stderr,
                    "\n\n *** clk_check_increments(): Bad increment!  Panicking!\n\n");
            exit(103);
         }
         if (ep->next)
         {
            if (ep->next->prev != ep)
            {
               fprintf(stderr,
                       "\n\n *** clk_check_increments(): Edges not connected properly!  Panicking!\n\n");
               exit(110);
            }
         }
         ep = ep->next;
      }
      if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN) && !equal_ivector(pos, start))
      {
         sub_ivector(pos, pos, knot->loffset);
         sub_ivector(start, start, knot->loffset);
         fprintf(stderr,
                 "\n\n *** clk_check_increments(): Component fails to close properly! (%d, %d, %d) != (%d, %d, %d) Panicking!\n\n",
                 start [0], start [1], start [2], pos [0], pos [1], pos [2]);

         clk_check_lattice(knot);
         exit(109);
      }
      comp = comp->next;
   }

}



#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

int clk_direction(ivector incr)
{
   // return the direction associated with an increment
   switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
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

bool clk_check_for_edge_hits(CubicLatticeKnotPtr knot, int dir, ivector test_location)
{
   // we should recentre the knot when an edge hit occurs
   // edge hits don't typically occur except for very small knots
   switch (dir)
   {
      case MOVE_NORTH:
         if (test_location [1] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_SOUTH:
         if (test_location [1] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_EAST:
         if (test_location [0] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_WEST:
         if (test_location [0] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_UP:
         if (test_location [2] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_DOWN:
         if (test_location [2] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
   }

   return false;
}

void clk_get_statistics(CubicLatticeKnotPtr clkp, int &n1, int &n2, int &n3, int &n4)
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

         else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
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

EdgePtr is_tight_clasp(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->dir != ep->prev->dir) return (EdgePtr) NULL;
   if (!perp [ep->dir][ep->next->dir]) return (EdgePtr) NULL;
   if (!anti [ep->prev->prev->dir][ep->next->dir]) return (EdgePtr) NULL;
   ivector loc;
   add_ivector(loc, ep->start, ep->next->increment);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   ivector incr;
   copy_ivector(incr, increment_NEWSUD [kross [ep->dir][ep->next->dir]]);

   ivector loc1;

   add_ivector(loc1, ep->start, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   sub_ivector(loc1, ep->start, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   add_ivector(loc1, loc, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   sub_ivector(loc1, loc, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep1 = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(loc, ep1->start))
         {
            if (ep1->dir != ep1->prev->dir) return (EdgePtr) NULL;
            if (!perp [ep1->dir][ep1->next->dir]) return (EdgePtr) NULL;
            if (!anti [ep1->prev->prev->dir][ep1->next->dir]) return (EdgePtr) NULL;
            ivector oloc;
            add_ivector(oloc, ep1->start, ep1->next->increment);
            if (!equal_ivector(oloc, ep->start)) return (EdgePtr) NULL;
            return ep1;
         }
         ep1 = ep1->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL; // should never get here
}

bool clk_check_path(char *path, CubicLatticeKnotPtr clkp, ivector loc, int value)
{
   ivector loco;
   int len = strlen(path);
   copy_ivector(loco, loc);
   for (int i = 0; i < len; i++)
   {
      add_ivector(loco, loco, increment_NEWSUD [(int) (path [i] - '0')]);
      if (clkp->lattice [lat(loco)] != value) return false;
   }
   return true;
}

EdgePtr is_parsite(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->dir != ep->prev->dir) return (EdgePtr) NULL;

   int xchange_dir;
   ivector loc;
   sub_ivector(loc, ep->start, clkp->loffset);

   switch (ep->dir)
   {
      case MOVE_NORTH:
      case MOVE_SOUTH:
         if (loc [0] > 0)
            xchange_dir = MOVE_EAST;
         else
            xchange_dir = MOVE_WEST;
         break;

      case MOVE_EAST:
      case MOVE_WEST:
         if (loc [1] > 0)
            xchange_dir = MOVE_NORTH;
         else
            xchange_dir = MOVE_SOUTH;
         break;

      default:
         return (EdgePtr) NULL;
   }


   char path [12];
   sprintf(path, "%d%d%d", xchange_dir, ep->dir, ep->dir);
   if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr) NULL;

   sprintf(path, "%d%d%d%d",
           MOVE_UP, ep->dir, xchange_dir, xchange_dir);
   if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr) NULL;

   ivector loc2;
   add_ivector(loc2, ep->start, increment_NEWSUD [xchange_dir]);
   add_ivector(loc2, loc2, increment_NEWSUD [xchange_dir]);
   sprintf(path, "%d%d%d%d",
           MOVE_DOWN, opposite [xchange_dir], opposite [xchange_dir], ep->dir);
   if (!clk_check_path(path, clkp, loc2, EMPTY)) return (EdgePtr) NULL;

   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep1 = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(loc2, ep1->start))
         {
            if (ep1->dir != ep1->prev->dir || ep1->dir != ep->dir) return (EdgePtr) NULL;
            //printf ("dir found is %s\n", clk_dir_name [xchange_dir]); fflush (stdout);
            return ep1;
         }
         ep1 = ep1->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL; // should never get here
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

void edge_info(char *b, EdgePtr ep)
{
   printf("%s, ID %d, loc %d %d %d, dir %d %d %d : %s, prev %d, next %d, ",
           b, ep->ID,
           ep->start [0] - 128, ep->start [1] - 128, ep->start [2] - 128,
           ep->increment [0], ep->increment [1], ep->increment [2],
           clk_dir_name [ep->dir],
           ep->prev->ID, ep->next->ID);
   if (ep->next->prev == ep)
      printf("Y");
   else
      printf("N");
   if (ep->prev->next == ep)
      printf("Y");
   else
      printf("N");

   ivector a;
   add_ivector(a, ep->start, ep->increment);
   if (equal_ivector(ep->next->start, a))
      printf("Y");
   else
      printf("N");
   add_ivector(a, ep->prev->start, ep->prev->increment);
   if (equal_ivector(ep->start, a))
      printf("Y");
   else
      printf("N");


   printf("\n");
   fflush(stdout);
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
   //  printf ("collapsing edge #%d\n", ep->ID);
   clkp->lattice [lat(ep->prev->start)] = EMPTY;
   clkp->lattice [lat(ep->next->start)] = EMPTY;

   copy_ivector(ep->prev->start, ep->prev->prev->start);
   add_ivector(ep->start, ep->prev->start, ep->prev->increment);

   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   delete_Edge(clkp, comp, ep->prev->prev);
   delete_Edge(clkp, comp, ep->next);
}

bool bfacf_strand_pass(CubicLatticeKnotPtr clkp)
{
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
   clkp->lattice [lat(ep1->start)] = OCCUPIED;
}

void panic_exit_parsite(char *c)
{
   fprintf(stderr, "\n\n *** panick exit!  %s\n", c);
   exit(343);
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, EdgePtr ep)
{
   int xchange_dir;
   EdgePtr epc;
   if (!(epc = is_parsite(clkp, ep))) return false;
   if (ep->dir != ep->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 1");
   if (ep->dir != epc->dir) panic_exit_parsite("bfacf_parsite_pass(): 2");
   if (epc->dir != epc->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 3");
   if (ep->dir == MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): 4");
   if (ep->dir == MOVE_DOWN) panic_exit_parsite("bfacf_parsite_pass(): 5");

   edge_info("ep", ep);

   ivector diff;
   bool found = false;
   sub_ivector(diff, epc->start, ep->start);
   for (xchange_dir = 0; xchange_dir < 6; xchange_dir++)
   {
      ivector incr;
      mult_ivector(incr, increment_NEWSUD [xchange_dir], 2);
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
   EdgePtr nep [16];
   int j = 2;
   for (int i = 0; i < 10; i++)
   {
      nep [j] = clkp->edgepool [clkp->nedges_total + i];
      nep [j]->comp = (void *) comp;
      j++;
      if (j == 10) j = 13;
   }
   nep [0] = epp;
   nep [1] = ep;
   nep [10] = epn;
   nep [11] = epcp;
   nep [12] = epc;
   nep [15] = epcn;

   int ep_orig_dir = ep->dir;
   // update positions 
   add_ivector(nep [1]->start, nep [0]->start, increment_NEWSUD [MOVE_UP]);
   add_ivector(nep [2]->start, nep [1]->start, increment_NEWSUD [ep->dir]);
   add_ivector(nep [3]->start, nep [2]->start, increment_NEWSUD [xchange_dir]);
   add_ivector(nep [4]->start, nep [3]->start, increment_NEWSUD [xchange_dir]);
   add_ivector(nep [5]->start, nep [4]->start, increment_NEWSUD [MOVE_DOWN]);
   add_ivector(nep [6]->start, nep [5]->start, increment_NEWSUD [MOVE_DOWN]);
   add_ivector(nep [7]->start, nep [6]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [8]->start, nep [7]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [9]->start, nep [8]->start, increment_NEWSUD [MOVE_UP]);

   add_ivector(nep [12]->start, nep [11]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [13]->start, nep [12]->start, increment_NEWSUD [ep->dir]);
   add_ivector(nep [14]->start, nep [13]->start, increment_NEWSUD [ep->dir]);

   // connect edges, adjust adjacency pointers and update lattice
   for (int i = 0; i < 15; i++)
   {
      if (i == 10) continue;
      join_edges(clkp, nep [i], nep [i + 1]);
   }

   // increase by 10 number of edges in comp and in knot
   comp->nedges += 10;
   clkp->nedges_total += 10;

   return true;
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp)
{
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
            ++kount;
            if (mark)
               mark(ep, epc);
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return kount;
}

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
         clkp->lattice [lat(ep->start)] = OCCUPIED;
         ep = ep->next;
      }
      comp = comp->next;
   }
}

//#define RGS_PARANOIA

bool recentre_knot_in_lattice(CubicLatticeKnotPtr knot)
{ // returns true if knot can't be recentered
   if (!knot) return true;

   ComponentCLKPtr comp;
   ivector min, max, mid;
   ivector diff;
   EdgePtr ep;

   clk_get_extent(min, max, knot);
   if (max [0] - min [0] > knot->max_range [0]) return true;
   if (max [1] - min [1] > knot->max_range [1]) return true;
   if (max [2] - min [2] > knot->max_range [2]) return true;

   // first clear the old path
   clear_lattice(knot);
   ep = knot->fcomp->first_edge;

#ifdef RGS_PARANOIA
   sub_ivector(diff, ep->start, knot->loffset);
   printf("before: %d %d %d   %d %d %d   %d %d %d\n",
           ep->start [0], ep->start [1], ep->start [2],
           knot->loffset [0], knot->loffset [1], knot->loffset [2],
           diff [0], diff [1], diff [2]);
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

   midpoint(mid, max, min);

   init_lattice(knot, mid);

#ifdef RGS_PARANOIA
   ep = knot->fcomp->first_edge;
   sub_ivector(diff, ep->start, knot->loffset);
   printf("after: %d %d %d   %d %d %d   %d %d %d\n\n",
           ep->start [0], ep->start [1], ep->start [2],
           knot->loffset [0], knot->loffset [1], knot->loffset [2],
           diff [0], diff [1], diff [2]);
   fflush(stdout);
#endif

   return false;
}

// currently, set_lattice(), fill_lattice(), invert_lattice(), dilate_lattice(), and clear_lattice()
// are used only by the KnotPlot driver to implement topological obstructions and cavities

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

   for (a [0] = offset_min [0]; a [0] <= offset_max [0]; a [0]++)
      for (a [1] = offset_min [1]; a [1] <= offset_max [1]; a [1]++)
         for (a [2] = offset_min [2]; a [2] <= offset_max [2]; a [2]++)
            clkp->lattice [lat(a)] = value;
}

void fill_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max)
{
   set_lattice(clkp, (char) OCCUPIED, min, max);
}

void invert_lattice(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return;
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (clkp->lattice [i] == OCCUPIED)
         clkp->lattice [i] = EMPTY;
      else
         clkp->lattice [i] = OCCUPIED;
   }
}

void dilate_lattice(CubicLatticeKnotPtr clkp)
{
   // expand occupied region to occupy neighboring cells
   if (!clkp) return;
   ivector pos;

   for (pos [2] = 1; pos [2] < LATTICE_SIZE - 1; pos [2]++)
   { // unnecessesary brace to make VC++ indenter happy
      for (pos [1] = 1; pos [1] < LATTICE_SIZE - 1; pos [1]++)
      {
         for (pos [0] = 1; pos [0] < LATTICE_SIZE - 1; pos [0]++)
         {
            if (clkp->lattice [lat(pos)] == OCCUPIED)
            {
               for (int i = 0; i < 6; i++)
               {
                  ivector pos2;
                  add_ivector(pos2, pos, increment_NEWSUD [i]);
                  if (clkp->lattice [lat(pos2)] == EMPTY)
                  {
                     clkp->lattice [lat(pos2)] = NEXTOCCUPIED;
                  }
               }
            }
         }
      }
   }
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (clkp->lattice [i] == NEXTOCCUPIED)
         clkp->lattice [i] = OCCUPIED;
   }
}

int bfacf_info_filled(CubicLatticeKnotPtr knot)
{
   int filled = 0;
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
      if (knot->lattice [i] == OCCUPIED)
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
         clkp->lattice [lat(ep->start)] = value;
         ep = ep->next;
      }
      if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
      {
         add_ivector(a, ep->start, ep->increment);
         clkp->lattice [lat(a)] = value;
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
   set_lattice(clkp, (char) EMPTY, min, max);
}

#define ABS(X)      ((X) < 0 ? -(X) : (X))

void bfacf_set_probabilities(ComponentCLKPtr comp, double pm2, double p0, double pp2)
{
   comp->p_minus2 = pm2;
   comp->p_0 = p0;
   comp->p_plus2 = pp2;

   // probably doesn't make much of a difference,
   // but precompute the following anyway:

   comp->p_4p2 = 4.0 * comp->p_plus2;
   comp->p_03p2 = comp->p_0 + 3.0 * comp->p_plus2;
   comp->p_m23p2 = comp->p_minus2 + 3.0 * comp->p_plus2;
   comp->p_2p0 = 2.0 * comp->p_0;
   comp->p_2p02p2 = 2.0 * comp->p_0 + 2.0 * comp->p_plus2;
}

void bfacf_set_probabilities(ComponentCLKPtr comp, double z)
{
   double p_plus2 = (z * z) / (1.0 + 3.0 * z * z);
   double p_0 = (1.0 + z * z) / (2.0 * (1.0 + 3.0 * z * z));
   double p_minus2 = 1.0 / (1.0 + 3.0 * z * z);
   comp->z = z;
   bfacf_set_probabilities(comp, p_minus2, p_0, p_plus2);
}

void bfacf_set_probabilities(CubicLatticeKnotPtr knot, double z)
{
   ComponentCLKPtr comp = knot->fcomp;
   while (comp)
   {
      bfacf_set_probabilities(comp, z);
      comp = comp->next;
   }
}

void bfacf_lattice_info(CubicLatticeKnotPtr knot)
{
   int min, max;
   int minloc, maxloc;
   minloc = maxloc = 0;
   min = max = knot->alt_lattice [0];
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (knot->alt_lattice [i] > max)
      {
         max = knot->alt_lattice [i];
         maxloc = i;
      }
      if (knot->alt_lattice [i] < min)
      {
         min = knot->alt_lattice [i];
         minloc = i;
      }
   }
   printf("%d, min value of %d at %d, max value of %d at %d\n", knot->nedges_total, min, minloc, max, maxloc);
}

void pexit(char *s)
{
   fprintf(stderr, "%s\n", s);
   fflush(stderr);
   exit(102);
}

bool clk_recombo_direct = true;

void clk_fix_incr(EdgePtr ep)
{
   sub_ivector(ep->increment, ep->next->start, ep->start);
   ep->dir = clk_direction(ep->increment);
   if (ep->dir == MOVE_INVALID)
      pexit("invalid increment detected");
}

bool perform_recombination_inverted(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2) {
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

    while (ep != ep2) {
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
    do {
        clk_fix_incr(ep);
        ep = ep->next;
        if (++kount > 16008)
            pexit("infinite loop B");
    } while (ep != ep2n);

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
   if (!anti [ep1->dir][ep2->dir])
      pexit("edges are not anti-parallel as expected!");


   ivector incr;
   sub_ivector(incr, ep1->next->start, ep2->start);
   if (!clk_check_increment(incr))
      pexit("invalid increment found!");

   ComponentCLKPtr comp1 = (ComponentCLKPtr) ep1->comp;
   ComponentCLKPtr comp2 = (ComponentCLKPtr) ep2->comp;

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
      comp2 = (ComponentCLKPtr) calloc(1, sizeof (ComponentCLK));
      comp2->minedges = clk_minedges;
      comp2->maxedges = clk_maxedges;
      memcpy(comp2, comp1, sizeof (ComponentCLK));
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
      comp2->next = (ComponentCLKPtr) NULL;
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
         ep->comp = (void *) comp;
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
      clkp->alt_lattice = (int *) calloc(LATTICE_TOTAL_SIZE, sizeof (int));
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
            clkp->alt_lattice [lat(ep1->start)] = ep1->locpool + 1;
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
               add_ivector(test_loc, ep1->start, increment_NEWSUD [turn [ep1->dir][dir]]);
               int IDp1 = clkp->alt_lattice [lat(test_loc)];
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
                     ep2 = clkp->edgepool [IDp1 - 1]->prev;
                     if (clk_recombo_limit && ep1->comp == ep2->comp && clkp->ncomps > 1)
                        ok = false;
                     if (ok && anti [ep1->dir][ep2->dir] && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2) >= clk_min_arc_length_distance &&
                             !(ep1->next == ep2->prev || ep1->prev == ep2->next))
                     {
                        ++kount;
                        mark(ep1, ep2);
                     }
                  }
                  else
                  {
                     ep2 = clkp->edgepool [IDp1 - 1];
                     if (ep1->dir == ep2->dir && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2) >= clk_min_arc_length_distance &&
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
            clkp->alt_lattice [lat(ep1->start)] = 0;
            ep1 = ep1->next;
         }
      }
      comp = comp->next;
   }

   return kount;
}
