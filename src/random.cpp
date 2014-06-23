/* 
 * File:   random.cpp
 * Author: kmo
 * 
 * Created on October 3, 2012, 3:40 PM
 */

#include "random.h"

#include <cmath>

pseudorandom::pseudorandom()
{
   current = 1;
   a = 16807;
   m = 2147483647;
   q = 127773;
   r = 2836;
   iset_gaussian = 0;
   gset_gaussian = 0.0;

   set_sRand_seed_to_clocktime();
}

pseudorandom::pseudorandom(int seed)
{
   current = 1;
   a = 16807;
   m = 2147483647;
   q = 127773;
   r = 2836;
   iset_gaussian = 0;
   gset_gaussian = 0.0;

   sRandSimple(seed);
}

pseudorandom::pseudorandom(const pseudorandom& orig) :
a(orig.a), current(orig.current), gset_gaussian(orig.gset_gaussian),
hi(orig.hi), iset_gaussian(orig.iset_gaussian), lo(orig.lo), m(orig.m),
q(orig.q), r(orig.r), test(orig.test) { }

pseudorandom::~pseudorandom() { }

#define TWOPI        6.28318530717958647693
#define ABS(X)      ((X) < 0 ? -(X) : (X))

int pseudorandom::RandSimple()
{
   hi = current / q;
   lo = current % q;
   test = a * lo - r * hi;
   current = (test > 0) ? test : test + m;
   return current;
}

/* Function rand_gaussian () returns a gaussian random variable */
/* with mean = 0 and sigma = 1. This routine is from Numerical  */

/* Recipes in C, page 217, function 'gasdev'                    */


double pseudorandom::rand_uniform()
{
//   return (double) RandSimple() / (double) 2147483647.0;
   return (double) RandSimple() / double(m);
}

double pseudorandom::rand_gaussian()
{
   double fac, r, v1, v2;

   if (iset_gaussian == 0)
   {
      do
      {
         v1 = 2.0 * rand_uniform() - 1.0;
         v2 = 2.0 * rand_uniform() - 1.0;
         r = v1 * v1 + v2 * v2;
      }
      while (r >= 1.0);
      fac = sqrt((double) (-2.0 * log(r) / r));
      gset_gaussian = v1 * fac;
      iset_gaussian = 1;
      return (double) (v2 * fac);
   }
   else
   {
      iset_gaussian = 0;
      return gset_gaussian;
   }
}

double pseudorandom::r_gaussian(double mean, double sigma)
{
   return mean + sigma * rand_gaussian();
}

///* 
//Function returns a random point uniformly distributed
//on a sphere.  Not really a great method to use.
//*/
//
//void rand_sphere (double *theta, double *phi) {
//  *theta = acos ((double) (2.0 * rand_uniform () - 1.0));
//  *phi   = TWOPI * rand_uniform ();
//}
//

/* function rand_integer (low, high) returns a random integer 
>= `low' and < `high' (NOTE: that's LESS THAN `high')
 */

int pseudorandom::rand_integer(int low, int high)
{
   return (int) ((double) (high - low) * rand_uniform() + (double) low);
}

double pseudorandom::rand_double(double low, double high)
{
   return ((high - low) * rand_uniform() + low);
}

//static int blurt = 0;

void pseudorandom::sRandSimple(int seed)
{
   int i;
   current = ABS(seed) + 2;

   iset_gaussian = 0;
   gset_gaussian = 0.0;

   /* 
   Cycle the random number generator a few times to get rid of
   any funny business.
    */

   for (i = 0; i < 51; i++) RandSimple();
}

void pseudorandom::sRandSimple(time_t t)
{
   sRandSimple((int) t);
}

void pseudorandom::set_sRand_seed_to_clocktime()
{
   time_t the_time;

   time(&the_time);

   sRandSimple(the_time);
}
