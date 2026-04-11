#include "legacy.h"

#include <cstdlib>
#include <cctype>
#include <cmath>
#include <ctime>

#define FALSE 0
#define TRUE 1

using namespace std;

#ifdef WIN32
#define ISNAN _isnan
#else
#define ISNAN std::isnan
#endif

// From random.cpp

#define TWOPI        6.28318530717958647693
#define ABS(X)      ((X) < 0 ? -(X) : (X))

//static int _rand_count = 0;
static int _Rand__SEED = 1,
        _Rand_a = 16807,
        _Rand_m = 2147483647,
        _Rand_q = 127773,
        _Rand_r = 2836,
        _Rand_lo, _Rand_hi, _Rand_test;

int RandSimple()
{
//   _rand_count++;
   _Rand_hi = _Rand__SEED / _Rand_q;
   _Rand_lo = _Rand__SEED % _Rand_q;
   _Rand_test = _Rand_a * _Rand_lo - _Rand_r * _Rand_hi;
   _Rand__SEED = (_Rand_test > 0) ? _Rand_test : _Rand_test + _Rand_m;
//   printf("%d %d ", _rand_count, _Rand__SEED);
   return _Rand__SEED;
}

/* Function rand_gaussian () returns a gaussian random variable */
/* with mean = 0 and sigma = 1. This routine is from Numerical  */
/* Recipes in C, page 217, function 'gasdev'                    */

static int iset_gaussian = 0;
static real gset_gaussian = 0.0;

real rand_uniform()
{
   return (real) RandSimple() / (real) 2147483647.0;
}

real rand_gaussian()
{
   real fac, r, v1, v2;

   if (iset_gaussian == 0)
   {
      do
      {
         v1 = 2.0 * rand_uniform() - 1.0;
         v2 = 2.0 * rand_uniform() - 1.0;
         r = v1 * v1 + v2 * v2;
      }
      while (r >= 1.0);
      fac = sqrt((real) (-2.0 * log(r) / r));
      gset_gaussian = v1 * fac;
      iset_gaussian = 1;
      return (real) (v2 * fac);
   }
   else
   {
      iset_gaussian = 0;
      return gset_gaussian;
   }
}

real r_gaussian(real mean, real sigma)
{
   return mean + sigma * rand_gaussian();
}

/* function rand_integer (low, high) returns a random integer
>= `low' and < `high' (NOTE: that's LESS THAN `high')
 */

int rand_integer(int low, int high)
{
   return (int) ((real) (high - low) * rand_uniform() + (real) low);
}

real rand_real(real low, real high)
{
   return ((high - low) * rand_uniform() + low);
}

static int blurt = 0;

void sRandSimple(int seed)
{
   int i;
   _Rand__SEED = ABS(seed) + 2;

   iset_gaussian = 0;
   gset_gaussian = 0.0;

   /*
   Cycle the random number generator a few times to get rid of
   any funny business.
    */

   for (i = 0; i < 51; i++) RandSimple();
}

void sRandSimple(time_t t)
{
   sRandSimple((int) t);
}

void set_sRand_seed_to_clocktime()
{
   time_t the_time;

   time(&the_time);

   sRandSimple(the_time);
}

// Reuben Brasher: the following added so that we can save the state of random
// number generator.

#include "pseudorandom.h"

pseudorandom saveRandomState()
{
   pseudorandom r;
   r.current = _Rand__SEED;
   r.a = _Rand_a;
   r.m = _Rand_m;
   r.q = _Rand_q;
   r.r = _Rand_r;
   r.lo = _Rand_lo;
   r.hi = _Rand_hi;
   r.test = _Rand_test;

   return r;
}

void copyRandomState(const pseudorandom& r)
{
   _Rand__SEED = r.current;
   _Rand_a = r.a;
   _Rand_m = r.m;
   _Rand_q = r.q;
   _Rand_r = r.r;
   _Rand_lo = r.lo;
   _Rand_hi = r.hi;
   _Rand_test = r.test;
}

// From sequence convert project vectory.cpp

vector3 iHat = {1.0, 0.0, 0.0};
vector3 jHat = {0.0, 1.0, 0.0};
vector3 kHat = {0.0, 0.0, 1.0};

vector3 xHat = {1.0, 0.0, 0.0};
vector3 yHat = {0.0, 1.0, 0.0};
vector3 zHat = {0.0, 0.0, 1.0};

vector3 ZERO_vect = {0.0, 0.0, 0.0};

#define EPSILON  0.00001

int round_to_long(real real_number)
{
   int sign = 1;
   if (real_number < 0.0)
   {
      sign = -1;
      real_number = -real_number;
   }
   return (sign * ((int) (real_number + 0.5)));
}

real dot(vector3 v, vector3 w)
{
   return (v [0] * w [0] + v [1] * w [1] + v [2] * w [2]);
}

real magnitude(vector3 v)
{
   return sqrt((v [0] * v [0] + v [1] * v [1] + v [2] * v [2]));
}

void cross(vector3 v, vector3 w, vector3 z)
{ // v = w X z
   v [0] = w [1] * z [2] - w [2] * z [1];
   v [1] = w [2] * z [0] - w [0] * z [2];
   v [2] = w [0] * z [1] - w [1] * z [0];
}

void sub_vector(vector3 v, vector3 w, vector3 z)
{
   v [0] = w [0] - z [0];
   v [1] = w [1] - z [1];
   v [2] = w [2] - z [2];
}

void normalize_vector(vector3 v, vector3 w)
{
   real length;

   if ((length = sqrt(dot(w, w))) != 0.0)
   {
      v [0] = w [0] / length;
      v [1] = w [1] / length;
      v [2] = w [2] / length;
   }
}

void set_vector(vector3 a, real x, real y, real z)
{
   a [0] = x;
   a [1] = y;
   a [2] = z;
}

void negate_ivector(ivector v, ivector w)
{
   v [0] = -w [0];
   v [1] = -w [1];
   v [2] = -w [2];
}

bool equal_ivector(ivector a, ivector b)
{
   if (a [0] != b [0]) return false;
   if (a [1] != b [1]) return false;
   if (a [2] != b [2]) return false;
   return true;
}

void sub_ivector(ivector v, ivector w, ivector z)
{
   v [0] = w [0] - z [0];
   v [1] = w [1] - z [1];
   v [2] = w [2] - z [2];
}

void add_ivector(ivector v, ivector w, ivector z)
{
   v [0] = w [0] + z [0];
   v [1] = w [1] + z [1];
   v [2] = w [2] + z [2];
}

void copy_ivector(ivector v, ivector w)
{
   v [0] = w [0];
   v [1] = w [1];
   v [2] = w [2];
}

void mult_ivector(ivector v, ivector w, double con)
{
   v [0] = round_to_long((double) w [0] * con);
   v [1] = round_to_long((double) w [1] * con);
   v [2] = round_to_long((double) w [2] * con);
}

void mult_ivector(ivector v, ivector w, int con)
{
   v [0] = w [0] * con;
   v [1] = w [1] * con;
   v [2] = w [2] * con;
}

void midpoint(ivector mid, ivector a, ivector b)
{
   add_ivector(mid, a, b);
   mult_ivector(mid, mid, 0.5);
}

void max_ivector(ivector v, ivector w, ivector z)
{
   v [0] = MAX(w [0], z [0]);
   v [1] = MAX(w [1], z [1]);
   v [2] = MAX(w [2], z [2]);
}

void min_ivector(ivector v, ivector w, ivector z)
{
   v [0] = MIN(w [0], z [0]);
   v [1] = MIN(w [1], z [1]);
   v [2] = MIN(w [2], z [2]);
}

void set_ivector(ivector a, int x, int y, int z)
{
   a [0] = x;
   a [1] = y;
   a [2] = z;
}


FILE *open_file(char *filename, char *access)
{
   FILE *fp = fopen(filename, access);
   if (!fp)
   {
      fprintf(stderr, " *** Can't open file `%s' for %s access!\n", filename, access);
      exit(1);
   }

   return fp;
}