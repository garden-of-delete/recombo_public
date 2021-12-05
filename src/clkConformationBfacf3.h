/* 
 * File:   clkConformationSeqconvert.h
 * Author: kmo
 *
 * Created on December 4, 2012, 3:08 PM
 */

#ifndef CLKCONFORMATIONSEQCONVERT_H
#define	CLKCONFORMATIONSEQCONVERT_H

#include <list>
#include <map>
#include <vector>

#include "clk.h"
#include "clkConformationAsList.h"
#include "clkConformationBfacf3.h"
#include "bfacfProbabilities.h"

#define DEFAULT_Z 0.20815
#define MAX_PRECOMPUTE_LENGTH 5000

//stores a set of BFACF probabilities
struct probs{
	double p_plus2;
	double p_minus2;
	double p_0;
};

class clkConformationBfacf3Component;

/**
 * Wrapper class for link implemented in bfacf3.
 */
class clkConformationBfacf3
{
public:
   /**
    * Constructor for single component link, a knot. The z-value will be set to
    * default value, and all probabilities will be calculated based on that
    * z-value. Sets random number generator with seed from system timer. Uses
    * default value for size of edgepool.
    * 
    * @param firstComponent a cubic lattice knot.
    */
   clkConformationBfacf3(const clk& firstComponent);
   
   /**
    * Constructor for two component link. The z-value will be set to
    * default value, and all probabilities will be calculated based on that
    * z-value. Sets random number generator with seed from system timer. Uses
    * default value for size of edgepool.
    * @param firstComponent a cubic lattice knot.
    * @param secondComponent a cubic lattice knot that should not intersect 
    * firstComponent at any vertex.
    */
   clkConformationBfacf3(const clk& firstComponent, const clk& secondComponent);
   
   /**
    * Destructor. Frees memory for edges, lattice and link and component 
    * structure.
    */
   virtual ~clkConformationBfacf3();

   /**
    * Not yet implemented.
    * @param component
    */
   void setKnot(const clk& component);
   
   /**
    * Not yet implemented.
    * @param component
    */
   void addComponent(const clk& component);

   /**
    * Number of components in link.
    * @return positive integer.
    */
   int size() const;

   /**
    * Populates a list with copies of components.
    * @param components will be cleared before being populated.
    */
   void getComponents(std::list<clkConformationAsList>& components) const;
   
   /**
    * Give access to one component of the link.
    * @param i index of component to access. Expected to be a non-negative integer.
    * @return component. Modifying the returned component will modify the
    * underlying bfacf3 component.
    */
   clkConformationBfacf3Component getComponent(int i);
   
   /**
    * Accesses z-value of component 0.
    * @return z-value.
    */
   double getZ() const;
   
   /**
    * Sets z-value and probabilities based on z-value for every component in 
    * link.
    * @param z should be a positive number smaller than critical z-value.
    */
   void setZ(double z);
   
   /**
    * Sets seed for random number generator. Warning, this will set seed for
    * all instances.
    * @param seed an integer.
    */
   void setSeed(int seed);
   
   /**
    * Attempt to perform a BFACF step on current conformation.
    */
   void step();
   
   /**
    * Attempt to perform a BFACF step on current conformation n times.
    */
   virtual void step(int n);
   
    /**
    * initializes an ordered hashmap of bfacf probability containing structures for a wide range of lengths
	*/
   void init_Q(double z, double q);

   /**
    * Attempt to perform a BFACF step on current conformation using Q
	*/
   void stepQ();

    /**
    * Attempt to perform a BFACF step on current conformation n times using Q.
    */
   //virtual void stepQ(int Q);

   virtual void stepQ(int q, double z);

   virtual void stepQ(long int c, int q, double z);

   /**
    * If conformation is a knot, a one-component link, then function will search
    * for all anti-parallel edges which have arcs between them of minimal length
    * at least minarclength and maximal arclength at most maxarclength.
    * 
    * </p >
    * 
    * If conformation is a two-component link, then function will search for
    * anti-parallel edges on opposite components as long as both components are
    * at least minarclength + 1 and at most maxarclength + 1.
    * 
    * </p >
    * 
    * In either of these two cases, calling this function will also prepare for
    * performing recombination between any of the recombination sites found.
    * 
    * @param minarclength must be at least five.
    * @param maxarclength must be at least minarclength.
    * @return non-negative integer, the count of pairs of edges found which
    * satisfy the minarclength and max arclength criterion. If return value is
    * zero, then recombination cannot be performed after calling this function.
    * If return value is a positive integer, then there will be that many
    * choices for of edge pairs for recombination.
    */

   int countRecomboSites(int minarclength, int maxarclength, int Sequence_type, int Recombo_type, int& Total_para_site, int& Total_anti_site, int& Para_site, int& Anti_site);
   int countRecomboSites(int minarclength, int maxarclength);
   
   /**
    * Perform a recombination between one of the edge pairs found by most 
    * recent call of countRecomboSites(int, int, char). Results will be unpredicable
    * if countRecomboSites(int, int, char) has not been called since the most recent
    * call of step() or step(int). Results will also be unpredicable if n is
    * greater than the return value of the most recent call to 
    * countRecomboSites(int, int, char).
    * @param n index of recombo site at which to perform recombination.
    */
   bool performRecombination(std::ostream& os, int Sequence_type, int Recombo_type, int component_num, int n = 0);
   bool performRecombination(int n = 0);







   /**
    * Will undo the most recent recombination operation performed by
    * performRecombination(int). Results will be unpredictable if no 
    * recombination has been performed, or if step() or step(int) have been 
    * called since the most recent call to performRecombination(int).
    */
   void undoRecombination(std::ostream& os, int Sequence_type, int Recombo_type, int component_num);
   void undoRecombination();
   probs* probMap;

    std::vector<threevector<int> > getChosenSite(int n);
    /**
     *given random choice from all site get the actual site from the list
     */

    bool dirIsParal(int site_choice);



protected:

   // Don't want default or copy constructors called, because they don't do anything now. Todo: make these constructors work
   clkConformationBfacf3();
   clkConformationBfacf3(const clkConformationBfacf3& orig);

private:
   class impl;
   impl *implementation;
   int q;
};

/**
 * Class giving access to bfacf probabilities stored in a bfacf3 component.
 */
class bfacfProbabilitiesFromComponent: public bfacfProbabilities
{
public:
   /**
    * Copies pointer to bfacf3 component.
    * @param orig
    */
   bfacfProbabilitiesFromComponent(const bfacfProbabilitiesFromComponent& orig);
   
   /**
    * Destructor.
    */
   virtual ~bfacfProbabilitiesFromComponent();

   virtual double p2() const;
   virtual double p0() const;
   virtual double m2() const;
   virtual double p4p2() const;
   virtual double p03p2() const;
   virtual double m23p2() const;
   virtual double p2p0() const;
   virtual double p2p02p2() const;
   
   friend class clkConformationBfacf3Component;
   
private:
   class impl;
   impl *implementation;
   
   // don't want to publicly call constructor
   bfacfProbabilitiesFromComponent(impl *implementation);

};

/**
 * Conformation whose size and vertices are components in bfacf3 implementation.
 * Wrapper class.
 */
class clkConformationBfacf3Component : public clk
{
public:
   virtual void getVertex(int index, threevector<int>& v) const;
   virtual void getVertices(std::list<threevector<int> >&) const;

   virtual int size() const;

   /**
    * Copy constructor copies point in original
    * @param orig
    */
   clkConformationBfacf3Component(const clkConformationBfacf3Component& orig);
   
   /**
    * Destructor does not free component owned by clkConformationBfacf3.
    */
   virtual ~clkConformationBfacf3Component();

   /**
    * Accesses z-value of bfacf3 component.
    * @return z-value.
    */
   double getZ() const;
   
   /**
    * Sets z-value and corresponding bfacf probabilities for one component only.
    * @param z should be a positive number less than the critical z-value.
    */
   void setZ(double z);
   
   /**
    * Accesses bfacf probabilities for this component.
    * @return 
    */
   bfacfProbabilitiesFromComponent getProbabilities() const;
   
protected:
   class impl;
   impl *implementation;

   // don't want constructor called, except by clkConformationBfacf3
   clkConformationBfacf3Component(impl *implementation);

   friend class clkConformationBfacf3;
};

//class cardinalDirection
//{
////   cardinalDirection(unsigned int);
////   cardinalDirection(unsigned short);
////   cardinalDirection(unsigned char);
//   cardinalDirection(const threevector<int>& v0, const threevector<int>& v1);
//   cardinalDirection(const threevector<int>& direction);
//   cardinalDirection(const cardinalDirection& orig);
//   virtual ~cardinalDirection();
//   
//   const std::string& name() const;
//   threevector<int> vector() const;
//   bool isopposite(const cardinalDirection& other) const;
//   
//private:
//   unsigned char impl;
//};

//class edgeBfacf3
//{
//   
//};
//
//class vertexBfacf3
//{
//   
//};
//
//class latticeBfacf3
//{
//   
//};

#endif	/* CLKCONFORMATIONSEQCONVERT_H */

