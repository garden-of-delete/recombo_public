/* 
 * File:   random.h
 * Author: kmo
 *
 * Created on October 3, 2012, 3:40 PM
 */

#ifndef RANDOM_H
#define	RANDOM_H

#include <ctime>

/**
 * A linear congruential pseudorandom number generator.
 */
class pseudorandom
{
public:
   /**
    * Default constructor sets seed using system timer.
    */
   pseudorandom();

   /**
    * Constructor initializes random number generator using seed.
    * @param seed an integer.
    */
   pseudorandom(int seed);

   /**
    * Copy constructor which duplicates state of random number generator. The
    * new random number generator will produce the same sequence as old one, 
    * starting in same place as old one.
    * @param orig
    */
   pseudorandom(const pseudorandom& orig);

   /**
    * Destructor.
    */
   virtual ~pseudorandom();

   /**
    * Uniform random real number.
    * @return uniformly between 0.0 and 1.0.
    */
   double rand_uniform();

   /**
    *  Generates a Gaussian random variable 
    * with mean = 0 and sigma = 1. This routine is from Numerical  
    * Recipes in C, page 217, function 'gasdev'
    * @return normal with mean 0.0 and standard deviation 1.0.
    */
   double rand_gaussian();

   /**
    * Generates a normal random variable with given mean and standard deviation.
    * @param mu a real number.
    * @param sigma a positive real number.
    * @return normal with mean mu and standard deviation sigma.
    */
   double r_gaussian(double mu, double sigma);

   /**
    * Uniformly generates a random integer within a given range.
    * @param low an integer.
    * @param high an integer greater than low.
    * @return an integer greater than or equal to low and less than or equal to
    * high.
    */
   int rand_integer(int low, int high);

   /**
    * Uniformly generates a random real within a given range.
    * @param low a real number.
    * @param high a real numerb greater than low.
    * @return a real number greater than or equal to low and less than or equal
    * to high.
    */
   double rand_double(double low, double high);

   /**
    * Initialize random number generator.
    * @param seed an integer.
    */
   void sRandSimple(int seed);


   /**
    * Initialize random number generator based on a time.
    * @param t will be converted to int.
    */
   void sRandSimple(time_t t);

   /**
    * Initialize random number generator using seed from system timer.
    */
   void set_sRand_seed_to_clocktime();

   /**
    * Default random number whose range depends on settings of linear 
    * congruential random number generator.
    * @return 
    */
   int RandSimple();

   /**
    * Save state of random number generator used by bfacf3 implementation.
    */
   friend pseudorandom saveRandomState();

   /**
    * Resets random number generator used by bfacf3 implementation to same state
    * as r.
    */
   friend void copyRandomState(const pseudorandom& r);

private:
   // variables used by simple random number generator
   int current,
   a,
   m,
   q,
   r,
   lo,
   hi,
   test;

   // variables used by Gaussian random number generator
   int iset_gaussian;
   double gset_gaussian;

   //   pseudorandom(bool);
};

pseudorandom saveRandomState();

#endif	/* RANDOM_H */

