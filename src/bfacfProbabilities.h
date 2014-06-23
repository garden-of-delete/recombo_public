/* 
 * File:   bfacfProbabilities.h
 * Author: kmo
 *
 * Created on January 19, 2013, 3:19 PM
 */

#ifndef BFACFPROBABILITIES_H
#define	BFACFPROBABILITIES_H

/**
 * Virtual class to access bfacf move probabilities and related quantities.
 */
class bfacfProbabilities
{
public:
   /**
    * Default constructor.
    */
   bfacfProbabilities();

   /**
    * Copy constructor.
    */
   bfacfProbabilities(const bfacfProbabilities& orig);

   /**
    * Destructor.
    */
   virtual ~bfacfProbabilities();

   /**
    * Plus two move.
    * @return Probability for accepting a +2 move.
    */
   virtual double p2() const = 0;

   /**
    * Zero move.
    * @return Probability for accepting a 0 move.
    */
   virtual double p0() const = 0;

   /**
    * Minus two move.
    * @return Probability for accepting a -2 move.
    */
   virtual double m2() const = 0;
   
   /**
    * Possibly precalculated.
    * @return 4 * p2()
    */
   virtual double p4p2() const = 0;
   
   /**
    * Possibly precalculated.
    * @return p0() + 3 * p2()
    */
   virtual double p03p2() const = 0;
   
   /**
    * Possibly precalculated.
    * @return m2() + 3 * p2()
    */
   virtual double m23p2() const = 0;
   
   /**
    * Possibly precalculated.
    * @return 2 * p0()
    */
   virtual double p2p0() const = 0;
   
   /**
    * Possibly precalculated.
    * @return 2 * p0() + 2 * p2()
    */
   virtual double p2p02p2() const = 0;

private:

};

#endif	/* BFACFPROBABILITIES_H */

