// Random number subroutines from KnotPlot

// Rob Scharein, April 2008

// Original algorithm from 
//  "Random number Generators: Good Ones Are Hard To Find",
//   Park, Stephen K., and Miller, Keith W.,
//   Communications of the ACM, 31 (10): 1192-1201., October, 1988.

// http://en.wikipedia.org/wiki/Park%E2%80%93Miller_RNG
// http://portal.acm.org/citation.cfm?id=63042

#include <math.h>
#include <time.h>
#include <stdlib.h>


#ifdef UNIX
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef WIN32
#include <process.h>
#endif

#include "random.h"

#define TWOPI        6.28318530717958647693
#define ABS(X)      ((X) < 0 ? -(X) : (X))

static int 
_Rand__SEED = 1,    
_Rand_a = 16807,
_Rand_m = 2147483647,
_Rand_q = 127773,
_Rand_r = 2836,
_Rand_lo, _Rand_hi, _Rand_test;

int RandSimple () {
  _Rand_hi = _Rand__SEED / _Rand_q;
  _Rand_lo = _Rand__SEED % _Rand_q;
  _Rand_test = _Rand_a * _Rand_lo - _Rand_r * _Rand_hi;
  _Rand__SEED = (_Rand_test > 0 ) ? _Rand_test : _Rand_test + _Rand_m;
  return _Rand__SEED;
}

/* Function rand_gaussian () returns a gaussian random variable */
/* with mean = 0 and sigma = 1. This routine is from Numerical  */
/* Recipes in C, page 217, function 'gasdev'                    */

static int iset_gaussian = 0;
static real gset_gaussian = 0.0;


real rand_uniform () {  
  return (real) RandSimple () / (real) 2147483647.0;
}

real rand_gaussian () {
  real fac, r, v1, v2;
  
  if (iset_gaussian == 0) {
    do {
      v1 = 2.0 * rand_uniform () - 1.0;
      v2 = 2.0 * rand_uniform () - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);
    fac = sqrt ((real) (-2.0 * log (r) / r));
    gset_gaussian = v1 * fac;
    iset_gaussian = 1;
    return (real) (v2 * fac);
  }
  else {
    iset_gaussian = 0;
    return gset_gaussian;
  }
}  


real r_gaussian (real mean, real sigma) {
  return mean + sigma * rand_gaussian ();
}



/* function rand_integer (low, high) returns a random integer 
>= `low' and < `high' (NOTE: that's LESS THAN `high')
*/

int rand_integer (int low, int high) {  
  return (int) ((real) (high - low) * rand_uniform () + (real) low);
}

real rand_real (real low, real high) {
  return ((high - low) * rand_uniform () + low);
}

static int blurt = 0;


void sRandSimple (int seed) {
  int i;
  _Rand__SEED = ABS (seed) + 2;
  
  iset_gaussian = 0;
  gset_gaussian = 0.0;
  
  /* 
  Cycle the random number generator a few times to get rid of
  any funny business.
  */
  
  for (i = 0; i < 51; i++) RandSimple ();
}

void sRandSimple (time_t t) {
  sRandSimple ((int) t);
}

void set_sRand_seed_to_clocktime () {
  time_t the_time;
  
  time (&the_time); 

  sRandSimple (the_time);
}
