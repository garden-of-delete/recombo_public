/* 
 * File:   clkConformationBfacf3.cpp
 * Author: kmo
 * 
 * Created on December 4, 2012, 3:08 PM
 */

#include "clkConformationBfacf3.h"
#include "legacyBfacf.h"
#include "threevector.h"   //diwen 08/22

#include <cstdlib>
#include <iostream>
#include <list>
#include <vector>
#include <utility>
#include "stdlib.h"
#include "stdio.h"

using namespace std;

#define ARC1(I, J)  ((J) - (I) - 1)
#define ARC2(I, J, N) ((N) + (I) - (J) - 1)
#define DEFAULT_RECOMBO_MAX 108
//maximum length to precompute BFACF transition probabilities when using the init_Q function


ostream& operator<<(ostream& out, ivector v)
{
   return out << v[0] << " " << v[1] << " " << v[2];
}

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

// implementation of bfacfProbabilitiesFromComponent

typedef const ComponentCLK *const_ComponentCLKPtr;

class bfacfProbabilitiesFromComponent::impl
{
public:
   const_ComponentCLKPtr component;

   impl(const_ComponentCLKPtr component) : component(component) { }

   ~impl() { }
};

bfacfProbabilitiesFromComponent::bfacfProbabilitiesFromComponent(bfacfProbabilitiesFromComponent::impl *implementation) : implementation(implementation) { }

bfacfProbabilitiesFromComponent::bfacfProbabilitiesFromComponent(const bfacfProbabilitiesFromComponent& orig) : implementation(orig.implementation) { }

bfacfProbabilitiesFromComponent::~bfacfProbabilitiesFromComponent()
{
   delete implementation;
}

double bfacfProbabilitiesFromComponent::p2() const
{
   return implementation->component->p_plus2;
}

double bfacfProbabilitiesFromComponent::p0() const
{
   return implementation->component->p_0;
}

double bfacfProbabilitiesFromComponent::m2() const
{
   return implementation->component->p_minus2;
}

double bfacfProbabilitiesFromComponent::p4p2() const
{
   return implementation->component->p_4p2;
}

double bfacfProbabilitiesFromComponent::p03p2() const
{
   return implementation->component->p_03p2;
}

double bfacfProbabilitiesFromComponent::m23p2() const
{
   return implementation->component->p_m23p2;
}

double bfacfProbabilitiesFromComponent::p2p0() const
{
   return implementation->component->p_2p0;
}

double bfacfProbabilitiesFromComponent::p2p02p2() const
{
   return implementation->component->p_2p02p2;
}



// implementation of clkConformationBfacf3Component

class clkConformationBfacf3Component::impl
{
public:
   ComponentCLKPtr component;

   impl(ComponentCLKPtr component) : component(component) { }

   ~impl() { }
};

void clkConformationBfacf3Component::getVertex(int index, threevector<int>& v) const
{
   index %= size();
   EdgePtr edge = implementation->component->first_edge;
   for (int i = 0; i < index; i++)
   {
      edge = edge->next;
   }
   //   v.set(edge->start[0], edge->start[1], edge->start[2]);
   v = edge->start;
}

void clkConformationBfacf3Component::getVertices(std::list<threevector<int> >& vertices) const
{
   vertices.clear();
   int nvertices = implementation->component->nedges;
   EdgePtr edge = implementation->component->first_edge;
   for (int i = 0; i < nvertices; i++)
   {
      vertices.push_back(threevector<int>(edge->start));
   }
}

int clkConformationBfacf3Component::size() const
{
   return implementation->component->nedges;
}

clkConformationBfacf3Component::clkConformationBfacf3Component(const clkConformationBfacf3Component & orig) :
implementation(new clkConformationBfacf3Component::impl(orig.implementation->component)) { }

clkConformationBfacf3Component::clkConformationBfacf3Component(clkConformationBfacf3Component::impl * implementation) :
implementation(implementation) { }

clkConformationBfacf3Component::~clkConformationBfacf3Component()
{
   // clkConformationBfacf3Component owns implementation, but none of the data in implementation
   delete implementation;
}

double clkConformationBfacf3Component::getZ() const
{
   return implementation->component->z;
}

void clkConformationBfacf3Component::setZ(double z)
{
   bfacf_set_probabilities(implementation->component, z);
}

bfacfProbabilitiesFromComponent clkConformationBfacf3Component::getProbabilities() const
{
   return bfacfProbabilitiesFromComponent(new bfacfProbabilitiesFromComponent::impl(implementation->component));
}

//class impl;
//impl *implementation;

// implementation of clkConformationBfacf3

static char directionString[][8] = {
   "north",
   "east",
   "west",
   "south",
   "up",
   "down",
   "unknown"
};

static int oppositeDir[6] = {
   MOVE_SOUTH,
   MOVE_WEST,
   MOVE_EAST,
   MOVE_NORTH,
   MOVE_DOWN,
   MOVE_UP
};

bool oppositeDirection(int d1, int d2)
{
   if (d1 < 0 || d1 >= 6) return false;
   return d2 == oppositeDir[d1];
}

class clkConformationBfacf3::impl
{
   CubicLatticeKnotPtr clkp;

   impl(const clk& firstComponent)
   {
      // TODO: might want to have the option of setting poolsize to something other than default
      int poolsize = DEFAULT_POOLSIZE;

      // need to keep track of bounds
      ivector min, max;
      clkp = (CubicLatticeKnotPtr) calloc(1, sizeof (CubicLatticeKnot));
      ComponentCLKPtr comp = createFirstComponent(firstComponent, clkp,
              min, max);

      // TODO: throw an exception if something goes wrong in creating a component
      //      if (!comp) clkp = NULL;

      ivector range;
      sub_ivector(range, max, min);
      if (range [0] > LATTICE_SIZE - 3 || range [1] > LATTICE_SIZE - 3 || range [2] > LATTICE_SIZE - 3)
      {
         cerr << " *** Knot has excessive extent in one or more dimensions" << endl;
         // TODO: should throw an exception rather than simply returning.
         return;
      }

      if (poolsize < (108 * clkp->nedges_total) / 100 + 108)
      {
         cerr << "poolsize too small" << endl;
         // TODO: should throw an exception rather than simply returning.
         //         free(clkp);
         //         return (CubicLatticeKnotPtr) NULL;
         return;
      }

      // ignore_mid to false
      bfacf_init_pool_and_lattice(clkp, false, poolsize, min, max);

      // set z to default
      bfacf_set_probabilities(clkp, DEFAULT_Z);

      // initialize machinery for recombo
      //      recombo_pair_index = 0;
      //      recombo_pair = (EdgePairPtr) calloc(DEFAULT_RECOMBO_MAX, sizeof (EdgePair));
      //      recombo_max = DEFAULT_RECOMBO_MAX;
   }

   impl(const clk& firstComponent, const clk& secondComponent)
   {
      // TODO: might want to have the option of setting poolsize to something other than default
      int poolsize = DEFAULT_POOLSIZE;

      // need to keep track of bounds
      ivector min, max;
      clkp = (CubicLatticeKnotPtr) calloc(1, sizeof (CubicLatticeKnot));
      ComponentCLKPtr comp1 = createFirstComponent(firstComponent, clkp,
              min, max);

      // TODO: throw an exception if something goes wrong in creating a component
      //      if (!comp1) clkp = NULL;

      ComponentCLKPtr comp2 = createNextComponent(secondComponent, clkp,
              min, max);

      // TODO: throw an exception if something goes wrong in creating a component
      //      if (!comp2) clkp = NULL;

      ivector range;
      sub_ivector(range, max, min);
      if (range [0] > LATTICE_SIZE - 3 || range [1] > LATTICE_SIZE - 3 || range [2] > LATTICE_SIZE - 3)
      {
         cerr << " *** Knot has excessive extent in one or more dimensions" << endl;
         // TODO: should throw an exception rather than simply returning.
         return;
      }

      if (poolsize < (108 * clkp->nedges_total) / 100 + 108)
      {
         cerr << "poolsize too small" << endl;
         // TODO: should throw an exception rather than simply returning.
         //         free(clkp);
         //         return (CubicLatticeKnotPtr) NULL;
         return;
      }

      // ignore_mid to false
      bfacf_init_pool_and_lattice(clkp, false, poolsize, min, max);

      // set z to default
      bfacf_set_probabilities(clkp, DEFAULT_Z);

      // initialize machinery for recombo
      //      recombo_pair_index = 0;
      //      recombo_pair = (EdgePairPtr) calloc(DEFAULT_RECOMBO_MAX, sizeof (EdgePair));
      //      recombo_max = DEFAULT_RECOMBO_MAX;
   }

   // TODO: right now this does nothing, so using this class creates a memory leak.

   ~impl()
   {
      //      cout << "Freeing link." << endl;

      // TODO: check to make sure that all resources are actually deallocated

      // free resources used for recombination
      //      if (recombo_pair) free(recombo_pair);

      if (!clkp)
      {
         // nothing to do
         return;
      }

      // TODO: Should we close the file here?
      //      if (clkp->fpout)
      //      {
      //         close(clkp->fpout);
      //         free(clkp->fpout);
      //      }

      // free the lattice
      if (clkp->lattice) free(clkp->lattice);
      // TODO: is alt_lattice actually used anywhere?
      if (clkp->alt_lattice) free(clkp->alt_lattice);

      // free the edgepool and the edges
      if (clkp->edgepool)
      {
         for (int i = 0; i < clkp->poolsize; i++)
         {
            if (clkp->edgepool[i]) free(clkp->edgepool[i]);
         }
         free(clkp->edgepool);
      }

      // free the components
      if (clkp->fcomp)
      {
         ComponentCLKPtr comp = clkp->fcomp;
         ComponentCLKPtr next;
         for (int i = 0; i < clkp->ncomps; i++)
         {
            next = comp->next;
            if (comp) free(comp);
            else break;
            comp = next;
         }
      }

      // finally free the link itself
      free(clkp);
   }

   // TODO: add implementations for next three functions

   void setKnot(const clk& component) { }

   void addComponent(const clk& component) { }

   void getComponents(std::list<clkConformationAsList>& components) const
   {
      // populate list with empty components
      components.clear();
      //      components.insert(components.begin(), clkp->ncomps, clkConformationAsList());

      ComponentCLKPtr component = clkp->fcomp;
      // iterate through the components
      for (int i = 0; i < clkp->ncomps; i++)
      {
         clkConformationBfacf3Component clkComponent(new clkConformationBfacf3Component::impl(component));
         components.push_back(clkComponent);
         component = component->next;
      }
   }

   int size() const
   {
      return clkp->ncomps;
   }

   ComponentCLKPtr createFirstComponent(const clk& nextComponent,
           CubicLatticeKnotPtr clkp,
           ivector min, ivector max)
   {
      ivector first_vertex, vertex, last_vertex;
      list<threevector<int> > vertices;
      nextComponent.getVertices(vertices);
      list<threevector<int> >::const_iterator i = vertices.begin();

      // remember the first vertex for later
      first_vertex[0] = last_vertex[0] = i->getX();
      first_vertex[1] = last_vertex[1] = i->getY();
      first_vertex[2] = last_vertex[2] = i->getZ();

      // create a new component
      ComponentCLKPtr comp = (ComponentCLKPtr) calloc(1, sizeof (ComponentCLK));
      comp->clkp = clkp;
      comp->minedges = clk_minedges;
      comp->maxedges = clk_maxedges;

      // set comp to be first component of link
      clkp->fcomp = clkp->lcomp = comp;
      clkp->ncomps = 1;


      // for now link is bounded by the only vertex
      copy_ivector(min, first_vertex);
      copy_ivector(max, first_vertex);

      i++;
      while (i != vertices.end())
      {
//         cout << *i << endl;
         vertex[0] = i->getX();
         vertex[1] = i->getY();
         vertex[2] = i->getZ();

         max_ivector(max, max, vertex);
         min_ivector(min, min, vertex);
         if (!add_edge_to_knot(clkp, comp, last_vertex, vertex))
         {
            cerr << "bad increment of ("
                    << comp->last_edge->increment [0] << ", "
                    << comp->last_edge->increment [1] << ", "
                    << comp->last_edge->increment [2] << ")" << endl;
            // TODO: should throw an exception rather than simply returning nothing        
            return (ComponentCLKPtr) NULL;
         }

         copy_ivector(last_vertex, vertex);

         i++;
      }

      // close off component
      add_edge_to_knot(clkp, comp, vertex, first_vertex);

      return comp;
   }

   ComponentCLKPtr createNextComponent(const clk& firstComponent,
           CubicLatticeKnotPtr clkp,
           ivector min, ivector max)
   {
      ivector first_vertex, vertex, last_vertex;
      list<threevector<int> > vertices;
      firstComponent.getVertices(vertices);
      list<threevector<int> >::const_iterator i = vertices.begin();

      // remember the first vertex for later
      first_vertex[0] = last_vertex[0] = i->getX();
      first_vertex[1] = last_vertex[1] = i->getY();
      first_vertex[2] = last_vertex[2] = i->getZ();

      // create a new component
      ComponentCLKPtr comp = (ComponentCLKPtr) calloc(1, sizeof (ComponentCLK));
      comp->clkp = clkp;
      comp->minedges = clk_minedges;
      comp->maxedges = clk_maxedges;

      // set comp to be next component of link
      clkp->lcomp->next = comp;
      comp->prev = clkp->lcomp;
      clkp->lcomp = comp;
      clkp->ncomps++;

      max_ivector(max, max, first_vertex);
      min_ivector(min, min, first_vertex);

      i++;
      while (i != vertices.end())
      {
         vertex[0] = i->getX();
         vertex[1] = i->getY();
         vertex[2] = i->getZ();

         max_ivector(max, max, vertex);
         min_ivector(min, min, vertex);
         if (!add_edge_to_knot(clkp, comp, last_vertex, vertex))
         {
            cerr << "bad increment of ("
                    << comp->last_edge->increment [0] << ", "
                    << comp->last_edge->increment [1] << ", "
                    << comp->last_edge->increment [2] << ")" << endl;
            // TODO: should throw an exception rather than simply returning nothing        
            return (ComponentCLKPtr) NULL;
         }

         copy_ivector(last_vertex, vertex);

         i++;
      }

      // close off component
      add_edge_to_knot(clkp, comp, vertex, first_vertex);

      return comp;

   }

   char* getDirectionString(int direction)
   {
      if (direction < 0 || direction >= 6) return directionString[6];
      return directionString[direction];
   }

   bool edgesUnitDistant(EdgePtr ep, EdgePtr eq)
   {
      // candidate edges may have opposite or same orientation
      if ((ep->dir == eq->dir) || oppositeDirection(ep->dir, eq->dir)) {
         //antiparallel distance
         ivector displacement_a;
         int dist_a;
         //parallel distance
         ivector displacement_p;
         int dist_p;
         // candidate edges must have beginning vertex distance one to ending vertex
         sub_ivector(displacement_a, ep->start, eq->next->start);
         sub_ivector(displacement_p, ep->start, eq->start);
         dist_a = abs(displacement_a[0]) + abs(displacement_a[1]) + abs(displacement_a[2]);
         dist_p = abs(displacement_p[0]) + abs(displacement_p[1]) + abs(displacement_p[2]);

         return (oppositeDirection(ep->dir, eq->dir) && dist_a == 1) || (ep->dir == eq->dir && dist_p == 1);
      }
      else
         return false;
   }

   bool edgesParallelAndUnitDistant(EdgePtr ep, EdgePtr eq)
   {
      // candidate edges must NOT have opposite orientation
      if (ep->dir!=eq->dir) return false;       
      
      ivector displacement;
      int dist;
      // candidate edges must have beginning vertex distance one to ending vertex
      sub_ivector(displacement, ep->start, eq->start);
      dist = abs(displacement[0]) + abs(displacement[1]) + abs(displacement[2]);
      //(debug) if (dist ==1) cout << ep->dir << eq->dir << dist;
      return dist == 1;
   }

   bool edgesAntiParallelAndUnitDistant(EdgePtr ep, EdgePtr eq)
   {
      // candidate edges must have opposite orientation
      if (!oppositeDirection(ep->dir, eq->dir)) return false;

      ivector displacement;
      int dist;
      // candidate edges must have beginning vertex distance one to ending vertex
      sub_ivector(displacement, ep->start, eq->next->start);
      dist = abs(displacement[0]) + abs(displacement[1]) + abs(displacement[2]);
      return dist == 1;
   }

   // making no attempt to optimize

   typedef pair<EdgePtr, EdgePtr> site;
   typedef std::vector<site> recomboSites;
   recomboSites sites;
   int lastRecombo;

    void searchForRecomboSitesBetweenComps() //CubicLatticeKnotPtr clkp) //, recomboSites& sites)
    {
       sites.clear();
       //      int count = 0;
       //      initRecombo();
       if (clkp->ncomps < 2)
       {
          // nothing more to do
          return;
       }

       ComponentCLK* c = clkp->fcomp;
       ComponentCLK* d = c->next;
       EdgePtr ep, eq;
       ep = c->first_edge;
       for (int i = 0; i < c->nedges; i++)
       {
          eq = d->first_edge;
          for (int j = 0; j < d->nedges; j++)
          {
             // candidate edges must have opposite orientation
             if (edgesAntiParallelAndUnitDistant(ep, eq))
             {
                //               count++;
                sites.push_back(site(ep, eq));
             }
             eq = eq->next;
          }
          ep = ep->next;
       }
       //      cout << "Number of sites recorded was " << count << endl;
    }

    void searchForRecomboSitesBetweenComps(int incoOrco, int stdOrvir, int vritualDir)//Used in mmc
   {
      sites.clear();
      //      int count = 0;
      //      initRecombo();
      if (clkp->ncomps < 2)
      {
         // nothing more to do
         return;
      }

      ComponentCLK* c = clkp->fcomp;
      ComponentCLK* d = c->next;
      EdgePtr ep, eq;
      ep = c->first_edge;
      for (int i = 0; i < c->nedges; i++)
      {
         eq = d->first_edge;
         for (int j = 0; j < d->nedges; j++)
         {
             if (stdOrvir == 2) { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,2,1), (1,2,0), (0,2,1), (0,2,0)
                 if (edgesUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
             }
             else if ((incoOrco == 1 && stdOrvir == 1) || (incoOrco == 0 && stdOrvir == 0)) { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,1), (0,0,1), (0,0,0)
                 if (edgesParallelAndUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
             }
             else { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,0,1), (1,0,0), (0,1)
                 if (edgesAntiParallelAndUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
             }

            eq = eq->next;
         }
         ep = ep->next;
      }
   }

   // making no attempt to optimize

   void searchForRecomboSitesOneComp(int minarclength, int maxarclength) //Only used in RunRecombo, has bugs
   {
      // for component length n, there will be two arcs between edges j > i
      // arc1 = j-i-1
      // arc2 = n+i-j-1
      // as j increases, arc1 increases and arc2 decreases
      //      int count = 0;
      //      initRecombo();
      //      this->clkp;
      sites.clear();
      int n = clkp->fcomp->nedges;
      // We must have 2 * minarclength <= n-2 <= 2 * maxarclength
      if (2 * minarclength > n - 2 || 2 * maxarclength < n - 2)
      {
         //printf("Length %d is not enough for minimum arclength %d and max %d.\n", n, minarclength, maxarclength);
         return;
      }

      ComponentCLK* c = clkp->fcomp;
      EdgePtr ep, eq;
      ep = c->first_edge;
      for (int i = 0; i < n; i++)
      {
         eq = ep->next;
         for (int j = i + 1; j < n; j++)
         {
            int arc1 = ARC1(i, j);
            int arc2 = ARC2(i, j, n);
            if (MIN(arc1, arc2) >= minarclength && MAX(arc1, arc2) <= maxarclength)
            {
               // candidate edges must have opposite orientation
               if (edgesAntiParallelAndUnitDistant(ep, eq))
               {
                  //                  count++;
                  sites.push_back(site(ep, eq));
               }

            }
            eq = eq->next;
         }
         ep = ep->next;
      }
      //      if (count > 0)
      //         cout << "Number of sites recorded was " << count << endl;
   }



   void searchForRecomboSitesOneComp(int minarclength, int maxarclength, int incoOrco, int stdOrvir, int vritualDir) //Used in mmc
   {
      sites.clear();
      int n = clkp->fcomp->nedges;

      // We must have 2 * minarclength <= n-2 <= 2 * maxarclength
      if (2 * minarclength > n - 2 || 2 * maxarclength < n - 2)
      {
         //printf("Length %d is not enough for minimum arclength %d and max %d.\n", n, minarclength, maxarclength);
         return;
      }

      ComponentCLK* c = clkp->fcomp;
      EdgePtr ep, eq;
      ep = c->first_edge;
      for (int i = 0; i < n; i++)
      {
         eq = ep->next;
         for (int j = i + 1; j < n; j++)
         {
            int arc1 = ARC1(i, j);
            int arc2 = ARC2(i, j, n);
            if (MIN(arc1, arc2) >= minarclength && MAX(arc1, arc2) <= maxarclength)
            {
               if (stdOrvir == 2) { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,2,1), (1,2,0), (0,2,1), (0,2,0)
                  if (edgesUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
               }
               else if ((incoOrco == 1 && stdOrvir == 1) || (incoOrco == 0 && stdOrvir == 0)) { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,1), (0,0,1), (0,0,0)
                  if (edgesParallelAndUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
               }
               else { //(inco_Or_co, std_Or_vir, virtual_Dir) = (1,0,1), (1,0,0), (0,1)
                  if (edgesAntiParallelAndUnitDistant(ep, eq))
                     sites.push_back(site(ep, eq));
               }
            }
            eq = eq->next;
         }
         ep = ep->next;
      }
   }

   friend class clkConformationBfacf3;
};


clkConformationBfacf3::clkConformationBfacf3(const clk & firstComponent) :
implementation(new clkConformationBfacf3::impl(firstComponent)) {
	probMap = new probs[MAX_PRECOMPUTE_LENGTH];
	q = 1;
}

clkConformationBfacf3::clkConformationBfacf3(const clk& firstComponent, const clk & secondComponent) :
implementation(new clkConformationBfacf3::impl(firstComponent, secondComponent)) {
	probMap = new probs[MAX_PRECOMPUTE_LENGTH];
	q = 1;
}	

clkConformationBfacf3::~clkConformationBfacf3()
{
   delete implementation;
   delete probMap;
}

void clkConformationBfacf3::setKnot(const clk & component)
{
   implementation->setKnot(component);
}

void clkConformationBfacf3::addComponent(const clk & component)
{
   implementation->addComponent(component);
}

void clkConformationBfacf3::getComponents(list<clkConformationAsList>& components) const
{
   implementation->getComponents(components);
}

clkConformationBfacf3Component clkConformationBfacf3::getComponent(int i)
{
   i %= implementation->clkp->ncomps;
   ComponentCLKPtr component = implementation->clkp->fcomp;
   for (int j = 0; j < i; j++)
   {
      component = component->next;
   }
   return clkConformationBfacf3Component(new clkConformationBfacf3Component::impl(component));
}

int clkConformationBfacf3::size() const
{
   return implementation->size();
}

double clkConformationBfacf3::getZ() const
{
   return implementation->clkp->fcomp->z;
}

void clkConformationBfacf3::setZ(double z)
{
	init_Q(z, q);
    bfacf_set_probabilities(implementation->clkp, z);
}

extern void sRandSimple(int seed);

void clkConformationBfacf3::setSeed(int seed)
{
   sRandSimple(seed);
}

void clkConformationBfacf3::step()
{
   perform_move(implementation->clkp);
}

void clkConformationBfacf3::step(int n)
{
   for (int i = 0; i < n; i++)
      perform_move(implementation->clkp);
}

void clkConformationBfacf3::init_Q(double z, double q) //these probabilities are suspect
{
	for (int i=4; i < MAX_PRECOMPUTE_LENGTH; i++){
		probMap[i].p_plus2 = (pow((i+2),(q-1))*(z * z)) / (pow(i,(q-1)) + 3.0*pow((i+2),q-1) * z * z);
		probMap[i].p_minus2 = pow(i,(q-1)) / (pow(i,(q-1)) + 3.0*pow((i+2),q-1) * z * z);
		probMap[i].p_0 = .5*(probMap[i].p_plus2 + probMap[i].p_minus2);
	}
}

void clkConformationBfacf3::stepQ(int current_q, double z)
{
	//perform_move_q(implementation->clkp); //will call the new step function(or possibly BE the new step function)
	if (current_q == 1){
		if (getZ() != z){
			bfacf_set_probabilities(implementation->clkp, z);
		}
		perform_move(implementation->clkp);
		return;
	}
	ComponentCLKPtr comp = implementation->clkp->fcomp;
	int n_comps = implementation->size();
	int n = 0;
	for (int i=0; i < n_comps; i++){
		n += getComponent(i).size();
	}
	/*
	//direct computation
	double p_plus2 = (pow((n+2),(q-1))*(z * z)) / (pow(n,(q-1)) + 3.0*pow((n+2),q-1) * z * z);
	double p_minus2 = pow(n,(q-1)) / (pow(n,(q-1)) + 3.0*pow((n+2),q-1) * z * z);
	double p_0 = .5*(p_plus2 + p_minus2);
	bfacf_set_probabilities(comp, p_minus2, p_0, p_plus2);*/
	
	if (getZ() != z || current_q != q){
		q = current_q;
		setZ(z);
	}
	if (n < MAX_PRECOMPUTE_LENGTH){
		bfacf_set_probabilities(comp, probMap[n].p_minus2, probMap[n].p_0, probMap[n].p_plus2); //changed from [4] to [n]
	}
	else{
		double p_plus2 = (pow((n+2),(q-1))*(z * z)) / (pow(n,(q-1)) + 3.0*pow((n+2),q-1) * z * z);
		double p_minus2 = pow(n,(q-1)) / (pow(n,(q-1)) + 3.0*pow((n+2),q-1) * z * z);
		double p_0 = .5*(p_plus2 + p_minus2);
		bfacf_set_probabilities(comp, p_minus2, p_0, p_plus2);
	}

	perform_move(implementation->clkp);
}

void clkConformationBfacf3::stepQ(long int c, int q, double z)
{
	for (long int i = 0; i < c; i++)
		stepQ(q, z);
		//perform_move_q(implementation->clkp); //see above
}


int clkConformationBfacf3::countRecomboSites(int minarclength, int maxarclength, int incoOrco, int stdOrvir, int virtualDir)
{
    //Sanity check
    if (maxarclength < minarclength)
    {
        cerr << "It does not make sense that maximum arc length should be less than minimum." << endl;
        exit(104);
    }
    // There must be at least three edges between edges of recombo site
    if (minarclength < 3)
    {
        cerr << "No possible recombo sites for minimum arc less than 3." << endl;
        exit(104);
    }
   //This version for use with orientation parameter
   if (size() == 1)
   {
      implementation->searchForRecomboSitesOneComp(minarclength, maxarclength, incoOrco, stdOrvir, virtualDir);
      return implementation->sites.size();
   }
   // here we assume that conformation is a two component link
   int sizeComp0, sizeComp1;
   sizeComp0 = getComponent(0).size();
   sizeComp1 = getComponent(1).size();
   if (sizeComp0 >= minarclength + 1 && sizeComp0 <= maxarclength + 1
          && sizeComp1 >= minarclength + 1 && sizeComp1 <= maxarclength + 1)
   {
      implementation->searchForRecomboSitesBetweenComps(incoOrco, stdOrvir, virtualDir); //minarclength, maxarclength);
      return implementation->sites.size();
   }
   return 0;
}

int clkConformationBfacf3::countRecomboSites(int minarclength, int maxarclength){ //only used by runRecombo.cpp
//Sanity check
    if (maxarclength < minarclength)
    {
        cerr << "It does not make sense that maximum arc length should be less than minimum." << endl;
        exit(104);
    }
    // There must be at least three edges between edges of recombo site
    if (minarclength < 3)
    {
        cerr << "No possible recombo sites for minimum arc less than 3." << endl;
        exit(104);
    }
    //This version for use with orientation parameter
    if (size() == 1)
    {
        implementation->searchForRecomboSitesOneComp(minarclength, maxarclength);
        return implementation->sites.size();
    }
    // here we assume that conformation is a two component link
    int sizeComp0, sizeComp1;
    sizeComp0 = getComponent(0).size();
    sizeComp1 = getComponent(1).size();
    if (sizeComp0 >= minarclength + 1 && sizeComp0 <= maxarclength + 1
        && sizeComp1 >= minarclength + 1 && sizeComp1 <= maxarclength + 1)
    {
        implementation->searchForRecomboSitesBetweenComps();
        return implementation->sites.size();
    }
    return 0;
}

class floatio //diwen 09/04
{
public:

    floatio() { }

    void write(int* arr, ostream& os) {
        float a = arr[0];
        float b = arr[1];
        float c = arr[2];

        write_float(&a, os);
        write_float(&b, os);
        write_float(&c, os);
    }

    void write(float* arr, ostream& os) {
        write_float(arr, os);
        write_float(arr + 1, os);
        write_float(arr + 2, os);
    }

    void write_int(int *val, ostream& os)
    {
        unsigned char a, b, c, d;

        a = (unsigned char) (0xff & (*val) >> 24);
        b = (unsigned char) (0xff & (*val) >> 16);
        c = (unsigned char) (0xff & (*val) >> 8);
        d = (unsigned char) (0xff & (*val));

        os.put(a);
        os.put(b);
        os.put(c);
        os.put(d);
    }

    void write_float(float *val, ostream& os){
        union {
            float f;
            int i;
        } x;
        x.f = *val;
        write_int (&x.i, os);
    }
};

void calculateNewCords(float *a, float *b, int *c, int *d , int *e, int *f, bool p){
    if (p) {
        a[0] = (c[0] + e[0]) * 0.5;
        a[1] = (c[1] + e[1]) * 0.5;
        a[2] = (c[2] + e[2]) * 0.5;
        b[0] = (d[0] + f[0]) * 0.5;
        b[1] = (d[1] + f[1]) * 0.5;
        b[2] = (d[2] + f[2]) * 0.5;
    }
    else{
        a[0] = (c[0] + f[0]) * 0.5;
        a[1] = (c[1] + f[1]) * 0.5;
        a[2] = (c[2] + f[2]) * 0.5;
        b[0] = (d[0] + e[0]) * 0.5;
        b[1] = (d[1] + e[1]) * 0.5;
        b[2] = (d[2] + e[2]) * 0.5;
    }
    if (c[0] == d[0] && d[0] == e[0] && e[0] == f[0]){
        a[0] += 0.5;
        b[0] -= 0.5;
    }
    else if (c[1] == d[1] && d[1] == e[1] && e[1] == f[1]){
        a[1] += 0.5;
        b[1] -= 0.5;
    }
    else {
        a[2] += 0.5;
        b[2] -= 0.5;
    }
}

bool clkConformationBfacf3::performRecombination(std::ostream& os, int incoOrco, int virtualDir, int component_num, int n)
{
   implementation->lastRecombo = n;
   clkConformationBfacf3::impl::site site = implementation->sites[n];
   //checks if the site contains parallel edges, need update ep2 after return from perform_recombination() to
   // make ep1 and ep2 parallel in the site
   // for the later use in undo_recombination()
   if (incoOrco == 1) //incoherent
   {
      if (site.first->dir == site.second->dir) //standard
      {
         EdgePtr newep2 = site.first->next;
         perform_recombination(implementation->clkp, site.first, site.second);
         implementation->sites[n].second = newep2;
         return true;
      }
      else // do virtual and write in binary
      {
         //binary input to os
          floatio floatio;
          float new_cordup[3], new_corddown[3];
          bool p = true;
          calculateNewCords(new_cordup, new_corddown, site.first->start, site.first->next->start, site.second->start, site.second->next->start, p);
          os.write("KnotPlot 1.0", 12);
          os.put((char)12);
          os.put((char)10);
          os.write("comp", 4);
          os.write("LOCF", 4);
          int Nvertices = implementation->clkp->nedges_total + 2;
          int num_bytes = 12 * Nvertices;
          floatio.write_int(&num_bytes, os);

          EdgePtr ep = site.first->next;
          while (ep != site.second){
              floatio.write(ep->start, os);
              ep = ep->next;
          }
          floatio.write(ep->start, os);
          if (virtualDir == 1)
              floatio.write(new_cordup, os);
          else
              floatio.write(new_corddown, os);
          floatio.write(site.first->start, os);
          ep = site.first->prev;
          while (ep != site.second){
              floatio.write(ep->start, os);
              ep = ep->prev;
          }
          if (virtualDir == 1)
              floatio.write(new_corddown, os);
          else
              floatio.write(new_cordup, os);
          os.write("endf", 4);
          os.put((char)10);
          return false;
      }
   }
   else  //coherent
   {
      if (site.first->dir != site.second->dir)
      {
         perform_recombination(implementation->clkp, site.first, site.second);
         return true;
      }
      else // do virtual
      {
         //binary input to os
         floatio floatio;
         float new_cordup[3], new_corddown[3];
         bool p = false;
         calculateNewCords(new_cordup, new_corddown, site.first->start, site.first->next->start, site.second->start, site.second->next->start, p);
         if (component_num == 1) {
            //writing the 1st comp
            os.write("KnotPlot 1.0", 12);
            os.put((char) 12);
            os.put((char) 10);
            os.write("comp", 4);
            os.write("LOCF", 4);
            int Nvertices1 = 0;
            EdgePtr ep1 = site.first;
            while (ep1 != site.second) {
               ep1 = ep1->next;
               Nvertices1++;
            }
            Nvertices1++;
            int num_bytes1 = 12 * Nvertices1;
            floatio.write_int(&num_bytes1, os);

            ep1 = site.first;
            while (ep1 != site.second) {
               ep1 = ep1->next;
               floatio.write(ep1->start, os);
            }
            //floatio.write(ep1->start, os);
            if (virtualDir == 1)
               floatio.write(new_cordup, os);
            else
               floatio.write(new_corddown, os);
            os.write("endf", 4);
            os.put((char) 10);
            //writing the 2nd comp
            os.write("KnotPlot 1.0", 12);
            os.put((char) 12);
            os.put((char) 10);
            os.write("comp", 4);
            os.write("LOCF", 4);
            int Nvertices2 = 0;
            EdgePtr ep2 = site.second;
            while (ep2 != site.first) {
               ep2 = ep2->next;
               Nvertices2++;
            }
            Nvertices2++;
            int num_bytes2 = 12 * Nvertices2;
            floatio.write_int(&num_bytes2, os);
            ep2 = site.second;
            while (ep2 != site.first) {
               ep2 = ep2->next;
               floatio.write(ep2->start, os);
            }
            if (virtualDir == 1)
               floatio.write(new_corddown, os);
            else
               floatio.write(new_cordup, os);
            os.write("endf", 4);
            os.put((char) 10);
            return false;
         }
         else //component_num = 2
         {
            os.write("KnotPlot 1.0", 12);
            os.put((char) 12);
            os.put((char) 10);
            os.write("comp", 4);
            os.write("LOCF", 4);
            int Nvertices = implementation->clkp->nedges_total + 2;
            int num_bytes = 12 * Nvertices;
            floatio.write_int(&num_bytes, os);
            EdgePtr ep = site.first;
            while (ep != site.first->prev){
                ep = ep->next;
                floatio.write(ep->start,os);
            }
            floatio.write(site.first->start, os);
            if (virtualDir == 1)
                floatio.write(new_cordup, os);
            else
                floatio.write(new_corddown, os);
            ep = site.second;
            while (ep != site.second->prev){
                ep = ep->next;
                floatio.write(ep->start, os);
            }
            floatio.write(site.second->start,os);
             if (virtualDir == 1)
                 floatio.write(new_corddown, os);
             else
                 floatio.write(new_cordup, os);
             os.write("endf", 4);
             os.put((char) 10);
             return false;
         }
      }
   }
}

bool clkConformationBfacf3::performRecombination(int n){   //only used by runRecombo.cpp
    implementation->lastRecombo = n;
    clkConformationBfacf3::impl::site site = implementation->sites[n];
    //checks if the site contains parallel edges, need update ep2 after return from perform_recombination() to
    // make ep1 and ep2 parallel in the site
    // for the later use in undo_recombination()
    if (site.first->dir == site.second->dir)
    {
        EdgePtr newep2 = site.first->next;
        perform_recombination(implementation->clkp, site.first, site.second);
        implementation->sites[n].second = newep2;
    }
    else
        perform_recombination(implementation->clkp, site.first, site.second);
}

std::vector<threevector<int> > clkConformationBfacf3::getChosenSite(int n){
   implementation->lastRecombo = n;
   clkConformationBfacf3::impl::site site = implementation->sites[n];
   std::vector<threevector<int> > recomboSites;


   recomboSites.push_back(threevector<int>(site.first->start));
   recomboSites.push_back(threevector<int>(site.first->next->start));
   recomboSites.push_back(threevector<int>(site.second->start));
   recomboSites.push_back(threevector<int>(site.second->next->start));
   std::string result = " ";
    for (int i = 0; i < 3; ++i) {
        result += std::to_string(site.first->start[i]);
        if(i != 2)
            result += ",";
    }
    result += " ";
    for (int i = 0; i < 3; ++i) {
        result += std::to_string(site.first->next->start[i]);
        if(i != 2)
            result += ",";
    }
    result += " ";
    for (int i = 0; i < 3; ++i) {
        result += std::to_string(site.second->start[i]);
        if(i != 2)
            result += ",";
    }
    result += " ";
    for (int i = 0; i < 3; ++i) {
        result += std::to_string(site.second->next->start[i]);
        if(i != 2)
            result += ",";
    }
    //result = site.first->start + " " + site.first->next->start + " " + site.second->start + " " + site.second->next->start;
    cout << result <<endl;
    return recomboSites;
}

bool clkConformationBfacf3::dirIsParal(int site_choice){
   if (implementation->sites[site_choice].first->dir == implementation->sites[site_choice].second->dir)
      return true;
   return false;
}

void clkConformationBfacf3::undoRecombination(std::ostream& os, int incoOrco, int virtualDir, int component_num)
{
   performRecombination(os, incoOrco, virtualDir, component_num, implementation->lastRecombo);
}

void clkConformationBfacf3::undoRecombination()  //only used by runRecombo.cpp
{
   performRecombination(implementation->lastRecombo);
}

// Don't want default or copy constructors called, because they don't do anything now. Todo: make these constructors work

clkConformationBfacf3::clkConformationBfacf3() { }

clkConformationBfacf3::clkConformationBfacf3(const clkConformationBfacf3 & orig) {
	
}

