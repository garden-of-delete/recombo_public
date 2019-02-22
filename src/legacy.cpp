#include "legacy.h"

#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cfloat>

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

static char comment [108];

vector3 iHat = {1.0, 0.0, 0.0};
vector3 jHat = {0.0, 1.0, 0.0};
vector3 kHat = {0.0, 0.0, 1.0};

vector3 xHat = {1.0, 0.0, 0.0};
vector3 yHat = {0.0, 1.0, 0.0};
vector3 zHat = {0.0, 0.0, 1.0};

vector3 ZERO_vect = {0.0, 0.0, 0.0};

#define EPSILON  0.00001

void output_it(char *s)
{
   // this function does something fancier within KnotPlot
   printf("%s", s);
}

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

void negate_vector(vector3 a, vector3 b)
{ // a = -b
   a [0] = -b [0];
   a [1] = -b [1];
   a [2] = -b [2];
}

int is_a_unit_vector(vector3 v)
{
   real mag;
   mag = magnitude(v);
   if (mag > 1.0 + EPSILON || mag < 1.0 - EPSILON)
      return FALSE;
   else
      return TRUE;
}

int is_a_zero_vector(vector3 v)
{
   if (v [0] > EPSILON) return FALSE;
   if (v [0] < -EPSILON) return FALSE;
   if (v [1] > EPSILON) return FALSE;
   if (v [1] < -EPSILON) return FALSE;
   if (v [2] > EPSILON) return FALSE;
   if (v [2] < -EPSILON) return FALSE;
   return TRUE;
}

int is_orthonormal_triad(vector3 u, vector3 v, vector3 w)
{
   // return TRUE iff (u, v, w) form an orthonormal triad
   vector3 x;

   if (!is_a_unit_vector(u)) return FALSE;
   if (!is_a_unit_vector(v)) return FALSE;
   if (!is_a_unit_vector(w)) return FALSE;
   cross(x, u, v);
   if (!is_a_unit_vector(x)) return FALSE;
   cross(x, v, w);
   if (!is_a_unit_vector(x)) return FALSE;
   cross(x, w, u);
   if (!is_a_unit_vector(x)) return FALSE;
   return TRUE;
}

void output_triad(vector3 u, vector3 v, vector3 w)
{
   vector3 uXv, vXw, wXu;

   sprintf(comment, "u = (%f, %f, %f) of length %f\n", u [0], u [1], u [2], magnitude(u));
   output_it(comment);
   sprintf(comment, "v = (%f, %f, %f) of length %f\n", v [0], v [1], v [2], magnitude(v));
   output_it(comment);
   sprintf(comment, "w = (%f, %f, %f) of length %f\n", w [0], w [1], w [2], magnitude(w));
   output_it(comment);

   cross(uXv, u, v);
   cross(vXw, v, w);
   cross(wXu, w, u);

   sprintf(comment, "mag (u X v) = %f, mag (v X w) = %f, mag (w X u) = %f\n",
           magnitude(uXv), magnitude(vXw), magnitude(wXu));
   output_it(comment);
}

void unit_vector(real theta, real phi, vector3 v)
{
   // creates a unit vector3 v from spherical coordinates theta, phi
   v [0] = sin(theta) * cos(phi);
   v [1] = sin(theta) * sin(phi);
   v [2] = cos(theta);
}

int equal_vector(vector3 a, vector3 b)
{
   if (a [0] != b [0]) return FALSE;
   if (a [1] != b [1]) return FALSE;
   if (a [2] != b [2]) return FALSE;
   return TRUE;
}

int equal_vector_within_epsilon(vector3 a, vector3 b)
{
   vector3 hypergumby;
   sub_vector(hypergumby, a, b);
   return is_a_zero_vector(hypergumby);
}

/*
Function returns a random point uniformly distributed
on a sphere.  Not really a great method to use.
 */

/*    // NOTE: use random_unit_vector () instead
void rand_sphere (real *theta, real *phi) {
 *theta = acos ((double) (2.0 * rand_uniform () - 1.0));
 *phi   = TWOPI * rand_uniform ();
}
 */

void random_unit_vector_cone(vector3 v, real alpha)
{ // full angle of the cone
   real theta, phi;

   //theta = acos (rand_real ((cos (alpha) + 1.0) * 0.5, 1.0));
   theta = acos(rand_real(cos(alpha * 0.5), 1.0));
   phi = TWOPI * rand_uniform();
   unit_vector(theta, phi, v);
}

void random_unit_vector(vector3 v)
{
   // creates a unit vector3 in a random direction
   // this method is fine in 3D, don't use in higher dimensions
   real dist;
   do
   { // loop until fall inside unit sphere
      v [0] = rand_real(-1.0, 1.0);
      v [1] = rand_real(-1.0, 1.0);
      v [2] = rand_real(-1.0, 1.0);
   }
   while (v [0] * v [0] + v [1] * v [1] + v [2] * v [2] > 1.0);
   normal_vector(v, v);
}

void random_vector(vector3 v, real mag)
{
   v [0] = rand_real(-mag, mag);
   v [1] = rand_real(-mag, mag);
   v [2] = rand_real(-mag, mag);
}

void zero_vector(vector3 v)
{
   v [0] = v [1] = v [2] = 0.0;
}

real dot(vector3 v, vector3 w)
{
   return (v [0] * w [0] + v [1] * w [1] + v [2] * w [2]);
}

real magnitude(vector3 v)
{
   return sqrt((v [0] * v [0] + v [1] * v [1] + v [2] * v [2]));
}

real magnitude_SQ(vector3 v)
{
   return (v [0] * v [0] + v [1] * v [1] + v [2] * v [2]);
}

void cross(vector3 v, vector3 w, vector3 z)
{ // v = w X z
   v [0] = w [1] * z [2] - w [2] * z [1];
   v [1] = w [2] * z [0] - w [0] * z [2];
   v [2] = w [0] * z [1] - w [1] * z [0];
}

void set_vector(vector3 a, real x, real y, real z)
{
   a [0] = x;
   a [1] = y;
   a [2] = z;
}

void set_fvector(fvector a, real x, real y, real z)
{
   a [0] = (freal) x;
   a [1] = (freal) y;
   a [2] = (freal) z;
}

void set_ivector4(ivector4 a, int x, int y, int z, int w)
{
   a [0] = x;
   a [1] = y;
   a [2] = z;
   a [3] = w;
}

void set_ivector2(ivector2 a, int x, int y)
{
   a [0] = x;
   a [1] = y;
}

void perp_vector(vector3 a, vector3 b)
{
   // creates a vector3 a perpendicular to b
   // because of `hairy ball theorem', this cannot be done
   // in a continuous fashion
   if (b [2] * b [2] > b [0] * b [0] + b [1] * b [1])
      cross(a, b, xHat);
   else
      cross(a, b, zHat);
}

void copy_vector(vector3 v, vector3 w)
{ // v = w
   v [0] = w [0];
   v [1] = w [1];
   v [2] = w [2];
}

void copy_ivector4(ivector4 v, ivector4 w)
{
   v [0] = w [0];
   v [1] = w [1];
   v [2] = w [2];
   v [3] = w [3];
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

void copy_from_fvector(vector3 v, fvector w)
{
   v [0] = (real) w [0];
   v [1] = (real) w [1];
   v [2] = (real) w [2];
}

void copy_to_fvector(fvector v, vector3 w)
{
   v [0] = (freal) w [0];
   v [1] = (freal) w [1];
   v [2] = (freal) w [2];
}

void add_vector(vector3 v, vector3 w, vector3 z)
{
   v [0] = w [0] + z [0];
   v [1] = w [1] + z [1];
   v [2] = w [2] + z [2];
}

void midpoint(vector3 mid, vector3 a, vector3 b)
{
   add_vector(mid, a, b);
   mult_vector(mid, mid, 0.5);
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

real length_segment(vector3 a, vector3 b)
{
   vector3 diff;
   sub_vector(diff, a, b);
   return magnitude(diff);
}

void rotate_vector_2D(vector3 a, vector3 b, real angle)
{
   real x, y;

   // rotate in 2D

   x = cos(angle) * b [0] - sin(angle) * b [1];
   y = sin(angle) * b [0] + cos(angle) * b [1];
   a [0] = x;
   a [1] = y;
}

void max_vector(vector3 v, vector3 w, vector3 z)
{
   v [0] = MAX(w [0], z [0]);
   v [1] = MAX(w [1], z [1]);
   v [2] = MAX(w [2], z [2]);
}

void min_vector(vector3 v, vector3 w, vector3 z)
{
   v [0] = MIN(w [0], z [0]);
   v [1] = MIN(w [1], z [1]);
   v [2] = MIN(w [2], z [2]);
}

void abs_vector(vector3 a, vector3 b)
{
   a [0] = ABS(b [0]);
   a [1] = ABS(b [1]);
   a [2] = ABS(b [2]);
}

void mult_vector(vector3 v, vector3 w, real con)
{
   v [0] = w [0] * con;
   v [1] = w [1] * con;
   v [2] = w [2] * con;
}

void interp_vector(vector3 a, vector3 b, vector3 c, real t)
{
   vector3 v;
   mult_vector(a, b, 1.0 - t);
   mult_vector(v, c, t);
   add_vector(a, a, v);
}

real interp_real(real b, real c, real t)
{
   return b * (1.0 - t) + t * c;
}

void transform_point(vector3 out, Matrix m, vector3 in)
{
   // out and in are column vectors, m is a matrix (usually a rotation matrix)
   // out = m in
#define M(row,col)  m[col*4+row]
   out [0] = M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3);
   out [1] = M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3);
   out [2] = M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3);
#undef M
}


void change_coord_system(vector3 new_a, vector3 a,
        vector3 xhat, vector3 yhat, vector3 zhat)
{
   vector3 dum;

   /* Two coordinate systems S and R,
      xhat yhat zhat are vectors expressed in S frame,
      xhat yhat zhat are orthonormal triad defining R frame,
      a is a vector3 expressed in the R frame,
      new_a is the same vector3 expressed in the S frame. */

   mult_vector(new_a, xhat, a [0]);
   mult_vector(dum, yhat, a [1]);
   add_vector(new_a, new_a, dum);
   mult_vector(dum, zhat, a [2]);
   add_vector(new_a, new_a, dum);
}

void change_coord_system_inv(vector3 new_a, vector3 a,
        vector3 xhat, vector3 yhat, vector3 zhat)
{

   /* Two coordinate systems S and R,
      xhat yhat zhat are vectors expressed in S frame,
      xhat yhat zhat are orthonormal triad defining R frame,
      a is a vector3 expressed in the S frame,
      new_a is the same vector3 expressed in the R frame. */

   new_a [0] = dot(a, xhat);
   new_a [1] = dot(a, yhat);
   new_a [2] = dot(a, zhat);
}

void new_normal_vector(vector3 v, vector3 w, int l, char *f)
{
   real length;

   if ((length = sqrt(dot(w, w))) != 0.0)
   {
      v [0] = w [0] / length;
      v [1] = w [1] / length;
      v [2] = w [2] / length;
   }
}

void make_orthonormal_triad(vector3 u, vector3 v, vector3 w)
{
   /* Makes an orthonormal right-handed triad given vector3 w. */

   perp_vector(u, w);
   new_normal_vector(u, u, __LINE__, __FILE__);
   cross(v, w, u);
}

void conv_to_xyz(vector3 vert, real radius, real theta, real phi)
{
   // vector3 `vert' is (radius, theta, phi) in spherical coordinates
   vert [0] = radius * sin(theta) * cos(phi);
   vert [1] = radius * sin(theta) * sin(phi);
   vert [2] = radius * cos(theta);
}

void random_orthonormal_triad(vector3 uhat, vector3 vhat, vector3 what)
{
   // creates a randomly oriented orthonormal triad
   vector3 u, v;
   real theta;
   random_unit_vector(what);
   make_orthonormal_triad(u, v, what);
   theta = rand_real(0.0, (real) TWOPI);
   mult_vector(u, u, cos(theta));
   mult_vector(v, v, sin(theta));
   add_vector(uhat, u, v);
   cross(vhat, what, uhat);
}

void shift_vector(vector3 v, int s)
{
   // component-wise shift of vector3
   vector3 a;
   copy_vector(a, v);
   v [0] = a [s % 3];
   v [1] = a [(s + 1) % 3];
   v [2] = a [(s + 2) % 3];
}

void compute_normal(vector3 n, vector3 v0, vector3 v1, vector3 v2)
{
   vector3 a, b;
   sub_vector(a, v2, v0);
   sub_vector(b, v1, v0);
   cross(n, b, a);
   new_normal_vector(n, n, __LINE__, __FILE__);
}

real diff_distance_SQ(vector3 a, vector3 b)
{
   // use this function if the magnitude of the difference is not needed
   // avoids taking an unnecessary square root
   vector3 diff;
   sub_vector(diff, a, b);
   return diff [0] * diff [0] + diff [1] * diff [1] + diff [2] * diff [2];
}

void perp_vector_wqed(vector3 a, vector3 b)
{
   // special function needed only in the animation done for WQED
   if (b [2] * b [2] > 20.0 * (b [0] * b [0] + b [1] * b [1]))
   {
      cross(a, b, xHat);
   }
   else
      cross(a, b, zHat);

   printf("xxxx");
   fflush(stdout);
}

void conv_cylindrical_to_xyz(vector3 vert,
        real radius, real theta, real z)
{
   vert [0] = radius * cos(theta);
   vert [1] = radius * sin(theta);
   vert [2] = z;
}

void conv_to_cylindrical(vector3 vert, real *radius, real *theta, real *z)
{
   *radius = sqrt(vert [0] * vert [0] +
           vert [1] * vert [1]);
   *theta = atan2(vert [1], vert [0]);
   *z = vert [2];
}

void conv_to_polar(vector3 vert, real *radius, real *theta, real *phi)
{
   *radius = sqrt(vert [0] * vert [0] +
           vert [1] * vert [1] +
           vert [2] * vert [2]);
   *phi = atan2(vert [1], vert [0]);
   if (*radius == 0.0)
      *theta = 0.0;
   else
      *theta = acos(vert [2] / *radius);
}

int is_right_of(vector3 p, vector3 base, vector3 dir)
{
   // assumes all vectors / points lie in xy-plane
   // given a ray with base point `base' and direction `dir',
   // returns TRUE iff point `p' is on the right hand side of the ray

   vector3 a, b;
   sub_vector(a, p, base);
   cross(b, a, dir);
   if (b [2] > 0.0)
      return TRUE;
   else
      return FALSE;
}

int inside_rect(vector3 a, vector3 b, real d, vector3 loc)
{
   // test point at location loc to see if it is
   // inside rectangle from a to b of width 2d

   // assumes vectors lie in xy-plane

   vector3 u; // axis of box
   vector3 v;
   vector3 v0, v1, v2, v3; // corners of box

   sub_vector(u, b, a);
   normal_vector(u, u);
   cross(v, zHat, u);
   mult_vector(v3, v, d);
   add_vector(v2, v3, b);
   add_vector(v3, v3, a);
   mult_vector(v0, v, -d);
   add_vector(v1, v0, b);
   add_vector(v0, v0, a);

   if (is_right_of(loc, v0, u)) return FALSE;
   if (is_right_of(loc, v1, v)) return FALSE;
   if (!is_right_of(loc, v3, u)) return FALSE;
   if (!is_right_of(loc, v0, v)) return FALSE;
   return TRUE;
}

void unit_vector_in_fan(vector3 f, vector3 g, point base, real angle)
{
   // create a new vector3 f which is perpendicular to g
   // if angle is zero, vector3 f is perpenicular to both g and the vector3 from the origin to the point base

   // first create a ugw coordinate system where g and (base - origin) lie in ug plane

   // base is a point on the unit sphere centered at the origin where the fan is anchored

   vector3 u, w;

   cross(w, base, g); // treat base as a vector3
   normal_vector(w, w);
   cross(u, g, w);

   mult_vector(u, u, cos(angle));
   mult_vector(w, w, sin(angle));
   add_vector(f, u, w);
}

void horizon_point(vector3 hp, vector3 cent, vector3 eye, vector3 up, real radius)
{
   // computes a horizon point

   // assumes `up' is a unit vector3

   real d, beta;
   vector3 a;

   sub_vector(a, cent, eye);
   d = magnitude(a);
   mult_vector(a, a, 1.0 / d); // make `a' a unit vector3
   beta = PI - acos(radius / d);
   mult_vector(hp, a, radius * cos(beta));
   mult_vector(a, up, radius * sin(beta));
   add_vector(hp, hp, a);
}

void get_sphere(vector3 centre, real &radius, vector3 *p)
{
   // given points p [0], p [1], p [2], and p [3] in general position,
   // compute the centre and radius of the sphere containing them

   Mat3x3 m;
   real d [3];

   for (int i = 0; i < 3; i++)
   {
      d [i] = 0.5 * (magnitude_SQ(p [i + 1]) - magnitude_SQ(p [0]));
      for (int j = 0; j < 3; j++)
      {
         m [i * 3 + j] = p [i + 1][j] - p [0][j];
      }
   }

   mat3X3inv(m);

   for (int k = 0; k < 3; k++)
   {
      centre [k] = dot(d, m + k * 3);
   }

   radius = sqrt(diff_distance_SQ(p [0], centre));
}

void mat3X3mult(Mat3x3 m1, Mat3x3 m2, Mat3x3 m3)
{
   int r, c, i;
   real sum;

   for (r = 0; r < 3; r++)
      for (c = 0; c < 3; c++)
      {
         sum = 0.0;
         for (i = 0; i < 3; i++)
            sum += (*(m2 + i + r * 3)) * (*(m3 + c + i * 3));
         *(m1 + c + r * 3) = sum;
      }
}

real mat3x3det(Mat3x3 mat)
{
   int r, c, r1, c1, r2, c2;
   real inv [9], det;

   for (r = 0; r < 3; r++)
   {
      for (c = 0; c < 3; c++)
      {
         r1 = (r + 1) % 3;
         c1 = (c + 1) % 3;
         r2 = (r + 2) % 3;
         c2 = (c + 2) % 3;

         inv [c + r * 3] = mat [c1 + r1 * 3] * mat [c2 + r2 * 3] -
                 mat [c1 + r2 * 3] * mat [c2 + r1 * 3];
      }
   }

   return mat [0] * inv [0] + mat [1] * inv [1] + mat [2] * inv [2];
}

void mat3X3inv(Mat3x3 mat)
{ /* Inefficient but suitable for where its needed */
   int r, c, r1, c1, r2, c2;
   real inv [9], det;

   for (r = 0; r < 3; r++)
   {
      for (c = 0; c < 3; c++)
      {
         r1 = (r + 1) % 3;
         c1 = (c + 1) % 3;
         r2 = (r + 2) % 3;
         c2 = (c + 2) % 3;

         inv [c + r * 3] = mat [c1 + r1 * 3] * mat [c2 + r2 * 3] -
                 mat [c1 + r2 * 3] * mat [c2 + r1 * 3];
      }
   }

   det = mat [0] * inv [0] + mat [1] * inv [1] + mat [2] * inv [2];

   if (det == 0.0)
   {
      fprintf(stderr, "Determinant is zero!\n");
      return;
   }

   for (r = 0; r < 3; r++)
      for (c = 0; c < 3; c++)
         mat [r + c * 3] = inv [c + r * 3] / det;

}

// From sequence convert project util.cpp


void init_increment(BitFilePtr bfp, FILE *fp);
bool write_increment(ivector incr, BitFilePtr bfp);
bool read_increment(ivector incr, BitFilePtr bfp);
int read_bit(BitFilePtr bfp);

int read_bit(BitFilePtr bfp);
void write_bit(int value, BitFilePtr bfp);

// KnotPlot 1.0 data files are stored using
// IRIX byte order ("big endian")

void read_unsigned_short(unsigned short *v, FILE *fp)
{
   unsigned char uca, ucb;
   fread(&ucb, 1, 1, fp);
   fread(&uca, 1, 1, fp);
   *v = (ucb << 8) | uca;
}

void write_int(int *val, FILE *fp)
{
   unsigned char a, b, c, d;

   a = (unsigned char) (0xff & (*val) >> 24);
   b = (unsigned char) (0xff & (*val) >> 16);
   c = (unsigned char) (0xff & (*val) >> 8);
   d = (unsigned char) (0xff & (*val));

   fwrite(&a, 1, 1, fp);
   fwrite(&b, 1, 1, fp);
   fwrite(&c, 1, 1, fp);
   fwrite(&d, 1, 1, fp);
}

void read_int(int *v, FILE *fp)
{
   unsigned char uc;
   *v = 0;
   for (int shift = 24; shift >= 0; shift -= 8)
   {
      fread(&uc, 1, 1, fp);
      *v |= uc << shift;
   }
}

void read_float(float *v, FILE *fp)
{

   union
   {
      int lg;
      float flt;
   } x;

   read_int(&x.lg, fp);
   *v = x.flt;
}

int needs_swapping()
{
   // returns TRUE on a "little endian" machine (Intel)
   // returns FALSE on a "big endian" machine (IRIX, Solaris)

   short super_gumby = 0x0001;
   char *swap_dat_thing = (char *) &super_gumby;
   return (int) *swap_dat_thing;
}

void read_double(double *v, FILE *fp)
{
   unsigned char *ucp;
   ucp = (unsigned char *) v;

   if (needs_swapping())
   {
      for (int i = 0; i < 8; i++)
         fread(ucp + 7 - i, 1, 1, fp);
   }
   else
   {
      for (int i = 0; i < 8; i++)
         fread(ucp++, 1, 1, fp);
   }
}

void read_fvector(fvector v, FILE *fp)
{
   read_float(v, fp);
   read_float(v + 1, fp);
   read_float(v + 2, fp);
}

void write_ivector(ivector v, FILE *fp)
{
   write_int(v, fp);
   write_int(v + 1, fp);
   write_int(v + 2, fp);
}

void read_ivector(ivector v, FILE *fp)
{
   read_int(v, fp);
   read_int(v + 1, fp);
   read_int(v + 2, fp);
}

/*
void sub_ivector (ivector v, ivector w, ivector z) {
  v [0] = w [0] - z [0];
  v [1] = w [1] - z [1];
  v [2] = w [2] - z [2];
}

void add_ivector (ivector v, ivector w, ivector z) {
  v [0] = w [0] + z [0];
  v [1] = w [1] + z [1];
  v [2] = w [2] + z [2];
}


void copy_ivector (ivector v, ivector w) {
  v [0] = w [0];
  v [1] = w [1];
  v [2] = w [2];
}

void mult_ivector (ivector v, ivector w, double con) {
  v [0] = round_to_long ((double) w [0] * con);
  v [1] = round_to_long ((double) w [1] * con);
  v [2] = round_to_long ((double) w [2] * con);
}

void mult_ivector (ivector v, ivector w, int con) {
  v [0] = w [0] * con;
  v [1] = w [1] * con;
  v [2] = w [2] * con;
}

void max_ivector (ivector v, ivector w, ivector z) {
  v [0] = MAX (w [0], z [0]);
  v [1] = MAX (w [1], z [1]);
  v [2] = MAX (w [2], z [2]);
}

void min_ivector (ivector v, ivector w, ivector z) {
  v [0] = MIN (w [0], z [0]);
  v [1] = MIN (w [1], z [1]);
  v [2] = MIN (w [2], z [2]);
}

void set_ivector (ivector a, int x, int y, int z) {
  a [0] = x;
  a [1] = y;
  a [2] = z;
}


int round_to_long (double real_number) {
  int sign = 1;
  if (real_number < 0.0) {
    sign = -1;
    real_number = -real_number;
  }
  return (sign * ((int) (real_number + 0.5)));
}

 */

void init_increment(BitFilePtr bfp, FILE *fp)
{
   bfp->ch = (char) 0;
   bfp->index = 0;
   bfp->fp = fp;
}

void write_bit(int b, BitFilePtr bfp)
{
   bfp->ch |= (b & 0x1) << bfp->index;
   if (++bfp->index == 8)
   {
      fwrite(&bfp->ch, 1, 1, bfp->fp);
      bfp->ch = 0;
      bfp->index = 0;
   }
}

int read_bit(BitFilePtr bfp)
{
   if (bfp->index % 8 == 0)
   {
      fread(&bfp->ch, 1, 1, bfp->fp);
      bfp->index = 0;
   }

   int value = bfp->ch >> bfp->index;
   bfp->index++;

   return 0x1 & value;
}

bool input_CUBE(int &num_vertices, FILE *fp, ivector *coord)
{
   char tag_name [5];

   if (fread(tag_name, 1, 4, fp) != 4) return false;

   tag_name [4] = (char) 0;
   if (strcmp(tag_name, "CUBE"))
   {
      fprintf(stderr, "bad tag >%s<\n", tag_name);
      fflush(stderr);
      return false;
   }
   int S;
   read_int(&S, fp);
   read_int(&num_vertices, fp);

   read_ivector(coord [0], fp);

   BitFile bf;
   init_increment(&bf, fp);
   ivector increment;

   for (int i = 1; i < num_vertices; i++)
   {
      if (!read_increment(increment, &bf)) return false;
      add_ivector(coord [i], coord [i - 1], increment);
   }
   if (!read_increment(increment, &bf)) return false; // ignore

   int numbits = num_vertices * 3;
   int numbytes = (numbits + 7) / 8;
   int num_pad_bits = numbytes * 8 - numbits;

   for (int pad = 0; pad < num_pad_bits; pad++)
   {
      if (read_bit(&bf) != 1)
      {
         fprintf(stderr, "bad bit seen!");
         fflush(stderr);
         return false;
      }
   }

   return true;
}

bool output_CUBE(int num_vertices, FILE *fp, ivector *coord)
{
   int numbits = num_vertices * 3;
   int numbytes = (numbits + 7) / 8;

   fprintf(fp, "CUBE");
   int S; // size of field is 16 bytes (number of vertices plus starting point) plus numbytes
   S = 16 + numbytes;
   write_int(&S, fp);


   // write number of vertices
   write_int(&num_vertices, fp);

   // write starting location

   ivector start;
   copy_ivector(start, coord [0]);
   write_ivector(start, fp);

   BitFile bf;
   init_increment(&bf, fp);

   for (int i = 1; i <= num_vertices; i++)
   {
      ivector increment;
      sub_ivector(increment, coord [i % num_vertices], coord [i - 1]);
      if (!write_increment(increment, &bf)) return false;
   }

   int num_pad_bits = numbytes * 8 - numbits;

   for (int pad = 0; pad < num_pad_bits; pad++)
      write_bit(1, &bf);

   return true;
}

#define MOVE_NORTH 0
#define MOVE_EAST  1
#define MOVE_WEST  2
#define MOVE_SOUTH 3
#define MOVE_UP    4
#define MOVE_DOWN  5

char dirname [6] = {'N', 'E', 'W', 'S', 'U', 'D'};

void write_increment(int value, BitFilePtr bfp)
{
   write_bit(value, bfp);
   write_bit(value >> 1, bfp);
   write_bit(value >> 2, bfp);
}

int read_increment(BitFilePtr bfp)
{
   int value = 0;
   for (int i = 0; i < 3; i++)
      value |= read_bit(bfp) << i;

   return value;
}

bool read_increment(ivector incr, BitFilePtr bfp)
{
   int dir;
   switch (dir = read_increment(bfp))
   {
      case MOVE_NORTH:
         set_ivector(incr, 0, 1, 0);
         break;
      case MOVE_EAST:
         set_ivector(incr, 1, 0, 0);
         break;
      case MOVE_WEST:
         set_ivector(incr, -1, 0, 0);
         break;
      case MOVE_SOUTH:
         set_ivector(incr, 0, -1, 0);
         break;
      case MOVE_DOWN:
         set_ivector(incr, 0, 0, -1);
         break;
      case MOVE_UP:
         set_ivector(incr, 0, 0, 1);
         break;
      default:
         fprintf(stderr, "bad move direction (%d) found!\n", dir);
         return false;
   }
   //  fprintf (stderr, "%c", dirname [dir]);

   return true;
}

bool write_increment(ivector incr, BitFilePtr bfp)
{
   int value;
   switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
   {
      case NORTH:
         value = MOVE_NORTH;
         break;
      case EAST:
         value = MOVE_EAST;
         break;
      case WEST:
         value = MOVE_WEST;
         break;
      case SOUTH:
         value = MOVE_SOUTH;
         break;
      case UP:
         value = MOVE_UP;
         break;
      case DOWN:
         value = MOVE_DOWN;
         break;
      default:
         fprintf(stderr, "bad increment (%d) (%d %d %d)\n",
                 (1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)),
                 incr [0], incr [1], incr [2]);
         return false;
   }
   write_increment(value, bfp);
   //  fprintf (stderr, "%c", dirname [value]);
   return true;
}

// From sequence convert project writhe.cpp

double writhe(double &ACN, int nvert, ivector *vert, double jitter){
    vector3 V1, V2, V3, V4;
    vector3 R13, R14, R24, R23, R34, R12;
    vector3 N1, N2, N3, N4;
    vector3 R34xR12;

    double gauss_xing_number_sum = 0.0;
    double space_writhe_sum = 0.0;
    double exact_acn_contrib;

    // first time thru the i loop requires special handling

    int last = nvert - 1;

    for (int i = 0; i < nvert; i++)
    {
        set_vector(V1,
                (double) vert [i][0] + rand_real(-jitter, jitter),
                (double) vert [i][1] + rand_real(-jitter, jitter),
                (double) vert [i][2] + rand_real(-jitter, jitter));

        set_vector(V2,
                (double) vert [(i + 1) % nvert][0] + rand_real(-jitter, jitter),
                (double) vert [(i + 1) % nvert][1] + rand_real(-jitter, jitter),
                (double) vert [(i + 1) % nvert][2] + rand_real(-jitter, jitter));


        sub_vector(R12, V2, V1);

        for (int j = i + 2; j < last; j++)
        {

            set_vector(V3,
                    (double) vert [j][0] + rand_real(-jitter, jitter),
                    (double) vert [j][1] + rand_real(-jitter, jitter),
                    (double) vert [j][2] + rand_real(-jitter, jitter));

            set_vector(V4,
                    (double) vert [(j + 1) % nvert][0] + rand_real(-jitter, jitter),
                    (double) vert [(j + 1) % nvert][1] + rand_real(-jitter, jitter),
                    (double) vert [(j + 1) % nvert][2] + rand_real(-jitter, jitter));

            sub_vector(R13, V3, V1);
            sub_vector(R14, V4, V1);
            sub_vector(R24, V4, V2);
            sub_vector(R23, V3, V2);
            sub_vector(R34, V4, V3);

            cross(N1, R13, R14);
            normalize_vector(N1, N1);
            cross(N2, R14, R24);
            normalize_vector(N2, N2);
            cross(N3, R24, R23);
            normalize_vector(N3, N3);
            cross(N4, R23, R13);
            normalize_vector(N4, N4);

            exact_acn_contrib = asin(dot(N1, N2)) + asin(dot(N2, N3)) +
                                asin(dot(N3, N4)) + asin(dot(N4, N1));

            if (ISNAN(exact_acn_contrib))
                continue;

            gauss_xing_number_sum += exact_acn_contrib;
            cross(R34xR12, R34, R12);


            if (dot(R34xR12, R13) > 0.0)
                space_writhe_sum += exact_acn_contrib;
            else
                space_writhe_sum -= exact_acn_contrib;


        }
        last = nvert;
    }

    double space_writhe = space_writhe_sum / TWOPI;
    ACN = gauss_xing_number_sum / TWOPI;
    return space_writhe;
}

double writhe_open(ivector* AB, ivector* CD) //Diwen 01/17/2019
{
    /*AB and CD are two straight lines in space with head and tail integer coordinates.*/
    double space_writhe = 0.0;
    vector3 A, B, C, D;
    vector3 R12, R13, R14, R23, R24, R34;
    vector3 N1, N2, N3, N4;
    vector3 R34xR12;

    set_vector(A, (double) AB[0][0], (double)AB[0][1], (double)AB[0][2]);
    set_vector(B, (double) AB[1][0], (double)AB[1][1], (double)AB[1][2]);
    set_vector(C, (double) CD[0][0], (double)CD[0][1], (double)CD[0][2]);
    set_vector(D, (double) CD[1][0], (double)CD[1][1], (double)CD[1][2]);

    sub_vector(R12, B, A);
    sub_vector(R13, C, A);
    sub_vector(R14, D, A);
    sub_vector(R23, C, B);
    sub_vector(R24, D, B);
    sub_vector(R34, D, C);

    cross(N1, R13, R14);
    cross(N2, R14, R24);
    cross(N3, R24, R23);
    cross(N4, R23, R13);

    normalize_vector(N1, N1);
    normalize_vector(N2, N2);
    normalize_vector(N3, N3);
    normalize_vector(N4, N4);

    double omega = asin(dot(N1, N2)) + asin(dot(N2, N3)) + asin(dot(N3, N4)) + asin(dot(N4, N1));
    cross(R34xR12, R34, R12);

    if (dot(R34xR12, R13) > 0.0)
        space_writhe += omega/TWOPI;
    else
        space_writhe -= omega/TWOPI;
    return space_writhe;
}

double radius_of_gyration(int nvert, ivector *vert)
{
   vector3 cofm, loc;
   zero_vector(cofm);
   for (int b = 0; b < nvert; b++)
   {
      set_vector(loc, (double) vert [b][0], (double) vert [b][1], (double) vert [b][2]);
      add_vector(cofm, cofm, loc);
   }
   mult_vector(cofm, cofm, 1.0 / (double) nvert);

   double sum = 0.0;
   for (int b = 0; b < nvert; b++)
   {
      set_vector(loc, (double) vert [b][0], (double) vert [b][1], (double) vert [b][2]);
      sum += diff_distance_SQ(loc, cofm);
   }

   return sqrt(sum / (double) nvert);
}

void bbox(int &volume, ivector range, int nvert, ivector *vert)
{
   ivector min, max;
   copy_ivector(min, vert [0]);
   copy_ivector(max, min);
   for (int b = 0; b < nvert; b++)
   {
      max_ivector(max, max, vert [b]);
      min_ivector(min, min, vert [b]);
   }
   sub_ivector(range, max, min);
   volume = range [0] * range [1] * range [2];
}

/*
if (isnan (exact_acn_contrib)) {
ivector i1, i2;
sub_ivector (i1, vert [(i + 1) % nvert], vert [i]);
sub_ivector (i2, vert [(j + 1) % nvert], vert [j]);

printf ("NaN: %d %d %d  %d %d %d  %d %d %d  %d %d %d\n",
vert [i][0], vert [i][1], vert [i][2],
vert [(i + 1) % nvert][0], vert [(i + 1) % nvert][1], vert [(i + 1) % nvert][2],
vert [j][0], vert [j][1], vert [j][2],
vert [(j + 1) % nvert][0], vert [(j + 1) % nvert][1], vert [(j + 1) % nvert][2]);
printf ("NaN: %d %d %d  %d %d %d\n",
i1 [0], i1 [1], i1 [2], i2 [0], i2 [1], i2 [2]);

}
 */

// From seqconvert.cpp main file

/*
seqconvert -h
for usage information

reads a knot sequence and writes sequences of various types

can operate as a pipe
see run-test.tcsh for examples of usage

 */

// available input / output formats
#define FORMAT_NOT_YET_SET 0
#define BINARY             1
#define PLAIN_TEXT         2
#define LENGTH             3
#define WRITHE             4
#define ENERGY             5
#define NEWSUD             6
#define ROG                7
#define BBOX               8
#define MULTIPLE_BINARY    9
#define PPM                10
#define ACN                11
#define UofS               12

// options
#define WRITE_NUM_VERTICES 1
#define CENTRE_KNOT        2
#define START_AT_ORIGIN    4
#define PROVIDE_INFO       8
#define NO_BLANK_LINE      16

#define OFFSET_SOMEHOW (CENTRE_KNOT | START_AT_ORIGIN)

int current_line = 1;
//char comment [108];
double writhe_mult = 1.0;
char *output_filename = (char *) NULL;
int energy_model;
double energy_parameter;
bool filter_energies = false;
double energy_scale = 1.0e21; // used only for values printed out, not used to compare energies

static ivector UofS_incr [7] = {
   {0, 0, 0}, // not used
   {1, 0, 0}, // 1   E
   {-1, 0, 0}, // 2   W
   {0, 1, 0}, // 3   N
   {0, -1, 0}, // 4   S
   {0, 0, 1}, // 5   U
   {0, 0, -1}
}; // 6   D

char north [4] = "3\n";
char east [4] = "1\n";
char west [4] = "2\n";
char south [4] = "4\n";
char up [4] = "5\n";
char down [4] = "6\n";

void put_a_knot_binary(FILE *fp, int flags);

void error_exit(char *s)
{
   fprintf(stderr, "\n  *** %s (line #%d)\n\n", s, current_line);
   exit(1);
}

void error_exit(char *s, char *t)
{
   fprintf(stderr, s, t);
   exit(1);
}

void error_exit(char *s, int t)
{
   fprintf(stderr, s, t);
   exit(1);
}

bool ends_in(char *ext, char *filename)
{
   // Test to see if string *t ends with extension *s, ignoring case

   if (!filename) return false;

   int sl, tl;
   sl = strlen(ext) - 1;
   tl = strlen(filename) - 1;
   if (sl > tl) return false;
   while (sl >= 0)
   {
      if (tolower((int) ext [sl]) != tolower((int) filename [tl])) return false;
      tl--;
      sl--;
   }
   return true;
}

#define MAXSIZE 30000

#ifdef HUGE
#undef HUGE
#endif

#define HUGE 2000000000

int min_length = 0; // knots must be of length >= min_length to get written to the output
int max_length = HUGE; // knots must be of length <= max_length to get written to the output

ivector global_max = {-HUGE, -HUGE, -HUGE};
ivector global_min = {HUGE, HUGE, HUGE};

ivector *coord;
ivector *coord_last;

bool save_last_knot = false;

int value [MAXSIZE * 3];
int num_vertices;
int num_converted = 0;

int num_rejected = 0, num_tested = 0;
bool debug_energy = false;

//bool reject_knot_energy_aux () {
//  if (!filter_energies) return false;
//
//  static double previous_energy;
//  static bool first_call = true;
//  static double X;
//
//  double current_energy = energy (coord, num_vertices);
//  num_tested++;
//
//  if (first_call) {
//    previous_energy = current_energy;
//    first_call = false;
//    X = 1.0 / (B * T);
//    return false;
//  }
//  if (debug_energy)
//    printf ("current %f, previous %f  ", current_energy * energy_scale, previous_energy * energy_scale);
//
//  if (current_energy > previous_energy) {
//    double p = exp ((previous_energy - current_energy) * X);  // probability of acceptance
//    if (debug_energy)
//      printf (" p = %f ", p );
//    if (rand_uniform () > p) {
//      if (debug_energy)
//        printf ("reject\n");
//      num_rejected++;
//      return true;                        // reject it
//    }
//  }
//  if (debug_energy)
//    printf ("accept\n");
//  previous_energy = current_energy;
//  return false;
//}
//
//bool reject_knot_energy () {
//  bool reject = reject_knot_energy_aux ();
//  if (!save_last_knot) return reject;
//  static int last_nvert;
//  if (reject) {
//    num_vertices = last_nvert;
//    for (int i = 0; i < num_vertices; i++)
//      copy_ivector (coord [i], coord_last [i]);
//  }
//  else {
//    last_nvert = num_vertices;
//    for (int i = 0; i < num_vertices; i++)
//      copy_ivector (coord_last [i], coord [i]);
//  }
//  return false;
//}

void offset_knot(ivector offset)
{
   if (offset [0] == 0 && offset [1] == 0 && offset [2] == 0) return;
   for (int i = 0; i < num_vertices; i++)
      sub_ivector(coord [i], coord [i], offset);
}

void process_knot(int flags)
{
   ivector max, min;
   copy_ivector(max, coord [0]);
   copy_ivector(min, coord [0]);

   for (int i = 0; i < num_vertices; i++)
   {
      max_ivector(max, max, coord [i]);
      min_ivector(min, min, coord [i]);
   }

   ivector offset = {0, 0, 0};

   switch (flags & OFFSET_SOMEHOW)
   {
      case CENTRE_KNOT:
         add_ivector(offset, min, max);
         mult_ivector(offset, offset, 0.5);
         break;
      case START_AT_ORIGIN:
         copy_ivector(offset, coord [0]);
         break;
   }

   offset_knot(offset);
   sub_ivector(max, max, offset);
   sub_ivector(min, min, offset);
   max_ivector(global_max, global_max, max);
   min_ivector(global_min, global_min, min);
}

/*
bool energy_allow (double &energy, int model, double param) {
switch (model) {
case 1:
return energy_model1 (energy, param, num_vertices, coord);
case 2:
return energy_model1 (energy, param, num_vertices, coord);
case 3:
return energy_model1 (energy, param, num_vertices, coord);
case 4:
return energy_model1 (energy, param, num_vertices, coord);
case 5:
return energy_model1 (energy, param, num_vertices, coord);
case 108:
return energy_model108 (energy, param, num_vertices, coord);
}
}
 */

void (*put_a_knot) (FILE *fp, int flags);

void put_a_knot_PPM(FILE *fp, int flags) { }

void put_a_knot_NEWSUD_or_UofS(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;
   process_knot(flags);
   if (north [0] == '3') // UofS format
      fprintf(fp, "%d %d %d\n", coord [0][0], coord [0][1], coord [0][2]);
   for (int i = 0; i < num_vertices; i++)
   {
      ivector incr;
      sub_ivector(incr, coord [(i + 1) % num_vertices], coord [i]);
      switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
      {
         case NORTH:
            fputs(north, fp);
            break;
         case EAST:
            fputs(east, fp);
            break;
         case WEST:
            fputs(west, fp);
            break;
         case SOUTH:
            fputs(south, fp);
            break;
         case UP:
            fputs(up, fp);
            break;
         case DOWN:
            fputs(down, fp);
            break;
      }
   }
   if (north [0] == '3') // UofS format
      fprintf(fp, "-111");
   fprintf(fp, "\n");
}

void put_a_knot_text(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   process_knot(flags);
   if (flags & WRITE_NUM_VERTICES)
      fprintf(fp, "%d ", num_vertices);
   for (int i = 0; i < num_vertices; i++)
      fprintf(fp, " %d %d %d", coord [i][0], coord [i][1], coord [i][2]);
   if (flags & NO_BLANK_LINE)
      fprintf(fp, "\n");
   else
      fprintf(fp, "\n\n");
   num_converted++;
}

/* added by Reuben Brasher
 For convenience, I wanted to list knot statistics with spaces between rather
 than new lines.
 */
char* separator = "\n";
//char* separator = " ";

void put_a_knot_length(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   process_knot(flags);

   fprintf(fp, "%d%s", num_vertices, separator);
   num_converted++;
}

int bbox_min_volume = HUGE;
int bbox_max_volume = 0;
ivector bbox_min_range = {HUGE, HUGE, HUGE};
ivector bbox_max_range = {0, 0, 0};

void put_a_knot_bbox(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   process_knot(flags);
   int volume;
   ivector range;

   bbox(volume, range, num_vertices, coord);
   bbox_min_volume = MIN(bbox_min_volume, volume);
   bbox_max_volume = MAX(bbox_max_volume, volume);

   for (int c = 0; c < 3; c++)
   { // component by component (x, y and z independently)
      bbox_min_range [c] = MIN(bbox_min_range [c], range [c]);
      bbox_max_range [c] = MAX(bbox_max_range [c], range [c]);
   }

   printf("%d %d %d %d\n", range [0], range [1], range [2], volume);

   num_converted++;
}

//int energy_number_skipped = 0;
//
//void put_a_knot_energy (FILE *fp, int flags) {
//  if (num_vertices > max_length || num_vertices < min_length) return;  // ignore this knot
////  if (reject_knot_energy ()) return;
//
//  process_knot (flags);
//
////  printf ("%f%s", energy_scale * energy (coord, num_vertices), separator);
//
//  num_converted++;
//}

#define DEFAULT_WRITHE_JITTER 0.0

double writhe_jitter = DEFAULT_WRITHE_JITTER;

void put_a_knot_writhe(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   double ignore;
   process_knot(flags);
   fprintf(fp, "%f%s", writhe_mult * writhe(ignore, num_vertices, coord, writhe_jitter), separator);
   num_converted++;
}

void put_a_knot_ACN(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   process_knot(flags);
   double acn;
   writhe(acn, num_vertices, coord, writhe_jitter);
   fprintf(fp, "%f%s", acn, separator);
   num_converted++;
}

void put_a_knot_rog(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   process_knot(flags);
   fprintf(fp, "%f%s", radius_of_gyration(num_vertices, coord), separator);
   num_converted++;
}

typedef struct
{
   FILE *fp;
   int kount;
} multFILE;

multFILE *mult_fp;

#define MAX_NUMB_SEQUENCES 10
int numb_sequences = 0;

void put_a_knot_multiple_binary(FILE *ignore, int flags)
{
   static int index = 0;
   static bool first_call = true;
   if (first_call)
   { // first call
      first_call = false;
      mult_fp = (multFILE *) calloc(MAX_NUMB_SEQUENCES, sizeof (multFILE));
      for (int i = 0; i < numb_sequences; i++)
      {
         sprintf(comment, "%s-%d.b", output_filename, i);
         mult_fp [i].fp = fopen(comment, "wb");
         if (!mult_fp [i].fp)
            error_exit("can't open file `%s' for write\n", comment);
      }
   }

   put_a_knot_binary(mult_fp [index].fp, flags);
   mult_fp [index].kount++;
   index = (index + 1) % numb_sequences;
}

void put_a_knot_binary(FILE *fp, int flags)
{
   if (num_vertices > max_length || num_vertices < min_length) return; // ignore this knot
   //  if (reject_knot_energy ()) return;

   if (num_vertices < 2)
      error_exit("absurd value for number of vertices!");


   process_knot(flags);

   if (!output_CUBE(num_vertices, fp, coord))
      error_exit("bad coordinate data!");

   num_converted++;

   current_line++;
}

bool (*get_a_knot) (FILE *fp);

bool get_a_knot_binary(FILE *fp)
{
   if (!input_CUBE(num_vertices, fp, coord)) return false;
   return true;
}

bool get_a_knot_text(FILE *fp)
{
   int index = 0;
   char ch;
   int sign = 1;
   int current_value = 0;
   bool digit_was_last_seen = false;
   int num_values;

   do
   {
      ch = fgetc(fp);
      switch (ch)
      {
         case EOF:
            return false;
         case ' ':
         case '\t':
         case '\n':
         case '\r':
            if (digit_was_last_seen)
            {
               value [index++] = current_value * sign;
               sign = 1;
               current_value = 0;
            }
            digit_was_last_seen = false;
            break;

         case '-':
            sign = -1;
            digit_was_last_seen = false;
            break;

         default:
            if (!isdigit((int) ch))
               error_exit("bad input seen!");
            current_value = 10 * current_value + ((int) (ch - '0'));
            digit_was_last_seen = true;
            break;
      }
   }
   while (ch != '\n' || index == 0);

   num_values = index;

   switch (num_values % 3)
   {
      case 0:
         num_vertices = num_values / 3;
         index = 0;
         break;
      case 1:
         num_vertices = (num_values - 1) / 3;
         if (num_vertices != value [0])
            error_exit("unexpected number of vertices!");
         index = 1;
         break;
      default:
         error_exit("wrong number of elements in input line!");
   }

   for (int i = 0; i < num_vertices; i++, index += 3)
   {
      set_ivector(coord [i], value [index], value [index + 1], value [index + 2]);
   }

   return true;
}

bool get_a_knot_UofS(FILE *fp)
{
   int x, y, z;
   fscanf(fp, "%d", &x);
   if (x == -999)
      return false;

   fscanf(fp, "%d%d", &y, &z);

   int index = 0;
   set_ivector(coord [index++], x, y, z);
   int value;
   do
   {
      fscanf(fp, "%d", &value);
      switch (value)
      {
         case -111:
            break;
         default:
            add_ivector(coord [index], coord [index - 1], UofS_incr [value]);
            index++;
      }
   }
   while (value != -111);

   // last coord [index] set is actually the same location as coord [0],
   // but coord [index] never gets used
   num_vertices = index - 1;
   return true;
}

bool get_a_knot_newsud(FILE *fp)
{
   int index = 0;
   char ch;

   set_ivector(coord [index++], 0, 0, 0);

   ivector incr = {NULL};

   do
   {
      ch = fgetc(fp);
      switch (ch)
      {
         case EOF:
            return false;
         case '\r':
         case '\n':
            break;
         case 'N':
            set_ivector(incr, 0, 1, 0);
            break;
         case 'E':
            set_ivector(incr, 1, 0, 0);
            break;
         case 'W':
            set_ivector(incr, -1, 0, 0);
            break;
         case 'S':
            set_ivector(incr, 0, -1, 0);
            break;
         case 'U':
            set_ivector(incr, 0, 0, 1);
            break;
         case 'D':
            set_ivector(incr, 0, 0, -1);
            break;
         default:
            return false;
      }
      if (isalpha((int) ch))
      {
         add_ivector(coord [index], coord [index - 1], incr);
         index++;
      }
   }
   while (ch != '\n' && index < MAXSIZE);

   num_vertices = index - 1;

   return true;
}


#define blurt if (wiki) fprintf (stderr, " ");

void usage(char *s)
{
   bool wiki = false;
   if (!strcmp(s, "wiki"))
      wiki = true;
   else if (strlen(s) > 2)
   {
      fprintf(stderr, "\n *** %s\n\n", s);
   }

   blurt;
   fprintf(stderr, "Sequence converter Version %s %s\n \n", __DATE__, __TIME__);
   blurt;
   fprintf(stderr, "Usage: seqconvert [options] inputfile outputfile\n");
   blurt;
   fprintf(stderr, "       or\n");
   blurt;
   fprintf(stderr, "       seqconvert [options] inputfile > outputfile\n");
   blurt;
   fprintf(stderr, "       or\n");
   blurt;
   fprintf(stderr, "       seqconvert -- [options] < inputfile > outputfile\n");
   blurt;
   fprintf(stderr, "       or\n");
   blurt;
   fprintf(stderr, "       some-sequence-producing-programme | seqconvert -- [options] > outputfile\n \n");
   blurt;
   fprintf(stderr, "       where the options are:\n");
   blurt;
   fprintf(stderr, "           -b            output in binary format\n");
   blurt;
   fprintf(stderr, "           -N            output in NEWSUD format\n");
   blurt;
   fprintf(stderr, "           -U            output in UofS format\n");
   blurt;
   fprintf(stderr, "           -t            output in plain text without number of vertices\n");
   blurt;
   fprintf(stderr, "           -v            output in plain text with number of vertices\n");
   blurt;
   fprintf(stderr, "           -n            don't write blank lines between knots\n");
   blurt;
   fprintf(stderr, "           -S PIECES     split a sequence file into PIECES sequences\n");
   blurt;
   fprintf(stderr, "           -V            same as -v and -n options combined\n");
   blurt;
   fprintf(stderr, "           -i            provide interesting information\n");
   blurt;
   fprintf(stderr, "           -q            run quietly\n");
   blurt;
   fprintf(stderr, "           -c            centre knot\n");
   blurt;
   fprintf(stderr, "           -s VALUE      set random number seed to VALUE (default is set to clock time)\n");
   blurt;
   fprintf(stderr, "           -o            translate knot so it starts at origin\n");
   blurt;
   fprintf(stderr, "           -r MIN MAX    specify a range for binning\n");
   blurt;
   fprintf(stderr, "           -f            filter knots according to energy\n");
   blurt;
   fprintf(stderr, "           -fr           filter knots according to energy and repeat knots\n");
   blurt;
   fprintf(stderr, "           -E            write only the energies to stdout\n");
   blurt;
   fprintf(stderr, "           -L            write only the lengths to stdout\n");
   blurt;
   fprintf(stderr, "           -W            write only the writhes to stdout\n");
   blurt;
   fprintf(stderr, "           -nW           write only the negative of the writhes to stdout\n");
   blurt;
   fprintf(stderr, "           -A            write only the average crossing number to stdout\n");
   blurt;
   fprintf(stderr, "           -R            write only the radii of gyration to stdout\n");
   blurt;
   fprintf(stderr, "           -j VALUE      set writhe jitter to VALUE (default is %f)\n", DEFAULT_WRITHE_JITTER);
   //  blurt;fprintf (stderr, "           -kappa KAPPA  set kappa value to KAPPA (default is %f)\n", DEFAULT_KAPPA);
   //  blurt;fprintf (stderr, "           -F  FACTOR    set f-factor value to FACTOR (default is %f)\n", DEFAULT_FUDGE_FACTOR);
   blurt;
   fprintf(stderr, "           -de           turn on energy debugging option (writes to stdout)\n");
   //  blurt;fprintf (stderr, "           -T TEMP       set temperature to TEMP degrees (default is %f)\n", DEFAULT_KELVINS);
   blurt;
   fprintf(stderr, "           -h            print this usage information\n");
   blurt;
   fprintf(stderr, "           --            read from standard input instead of trying to open a file\n \n");
   blurt;
   fprintf(stderr, "For more information see http://ewok.sfsu.edu/wiki/index.php/Sequence_Convert\n");
   fprintf(stderr, "\n");
   exit(102);
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

int legacy_main(int argc, char **argv)
{
   FILE *fpin;
   FILE *fpout = stdout;

   int flags = 0;

   bool read_from_stdin = false;
   int input_format = FORMAT_NOT_YET_SET;
   int output_format = FORMAT_NOT_YET_SET;

   int index = 1;

   bool write_only_energies = false;
   set_sRand_seed_to_clocktime();

   coord = (ivector *) calloc(MAXSIZE, sizeof (ivector));

   while (index < argc)
   {
      if (argv [index][0] != '-') break;

      if (!strcmp(argv [index], "-v"))
      {
         flags |= WRITE_NUM_VERTICES;
         output_format = PLAIN_TEXT;
      }

      else if (!strcmp(argv [index], "-V"))
      {
         flags |= (WRITE_NUM_VERTICES | NO_BLANK_LINE);
         output_format = PLAIN_TEXT;
      }

      else if (!strcmp(argv [index], "-n"))
      {
         flags |= NO_BLANK_LINE;
      }

      else if (!strcmp(argv [index], "-b"))
         output_format = BINARY;

      else if (!strcmp(argv [index], "-N"))
      {
         output_format = NEWSUD;
         sprintf(north, "N");
         sprintf(east, "E");
         sprintf(west, "W");
         sprintf(south, "S");
         sprintf(up, "U");
         sprintf(down, "D");
      }

      else if (!strcmp(argv [index], "-U"))
      {
         output_format = UofS;
      }

      else if (!strcmp(argv [index], "--"))
         read_from_stdin = true;

      else if (!strcmp(argv [index], "-t"))
         output_format = PLAIN_TEXT;

      else if (!strcmp(argv [index], "-c"))
         flags |= CENTRE_KNOT;

      else if (!strcmp(argv [index], "-p"))
         fprintf(stderr, "-p option no longer needed\n");

      else if (!strcmp(argv [index], "-o"))
         flags |= START_AT_ORIGIN;

      else if (!strcmp(argv [index], "-i"))
         flags |= PROVIDE_INFO;

      else if (!strcmp(argv [index], "-q"))
         flags &= ~PROVIDE_INFO;

      else if (!strcmp(argv [index], "-E"))
      {
         error_exit("option `%s' currently broken\n", argv [index]);
         energy_model = 108;
         output_format = ENERGY;
      }

      else if (!strcmp(argv [index], "-L"))
         output_format = LENGTH;

      else if (!strcmp(argv [index], "-A"))
         output_format = ACN;

      else if (!strcmp(argv [index], "-W"))
         output_format = WRITHE;

      else if (!strcmp(argv [index], "-nW"))
      {
         output_format = WRITHE;
         writhe_mult = -1.0;
      }

      else if (!strcmp(argv [index], "-R"))
         output_format = ROG;

      else if (!strcmp(argv [index], "-bb"))
         output_format = BBOX;

      else if (!strcmp(argv [index], "-ppm"))
         output_format = PPM;

      else if (!strcmp(argv [index], "-wiki"))
         usage("wiki");

      else if (!strcmp(argv [index], "-s"))
      {
         int seed = atoi(argv [index + 1]);
         fprintf(stderr, "Setting random number generator seed to %d.\n", seed);
         sRandSimple(seed);
         index++;
      }

      else if (!strcmp(argv [index], "-S"))
      {
         numb_sequences = atoi(argv [index + 1]);
         if (numb_sequences < 1 || numb_sequences > MAX_NUMB_SEQUENCES - 1)
            error_exit("absurd value of `%d' for number of output sequences\n", numb_sequences);
         output_format = MULTIPLE_BINARY;
         index++;
      }

      else if (!strcmp(argv [index], "-r"))
      {
         if (index + 2 >= argc)
            usage("needs two integers for length range");
         min_length = atoi(argv [index + 1]);
         max_length = atoi(argv [index + 2]);
         index += 2;
      }

      else if (!strcmp(argv [index], "-j"))
      {
         if (index + 1 >= argc)
            usage("needs a number for the writhe jitter");
         writhe_jitter = atof(argv [++index]);
      }

         //    else if (!strcmp (argv [index], "-kappa")) {
         //      if (index + 1 >= argc)
         //        usage ("expecting a value for kappa");
         //      energy_set_kappa (atof (argv [++index]));
         //    }
         //
         //    else if (!strcmp (argv [index], "-T")) {
         //      if (index + 1 >= argc)
         //        usage ("expecting a value for T");
         //      energy_set_T (atof (argv [++index]));
         //    }
         //
         //    else if (!strcmp (argv [index], "-F")) {
         //      if (index + 1 >= argc)
         //        usage ("expecting a value for the F factor");
         //      energy_set_fudge_factor (atof (argv [++index]));
         //    }

      else if (!strcmp(argv [index], "-f"))
      {
         error_exit("option `%s' currently broken\n", argv [index]);
         filter_energies = true;
      }

      else if (!strcmp(argv [index], "-fr"))
      {
         error_exit("option `%s' currently broken\n", argv [index]);
         filter_energies = true;
         save_last_knot = true;
      }

      else if (!strcmp(argv [index], "-h"))
         usage("");

      else if (!strcmp(argv [index], "-de"))
      {
         error_exit("option `%s' currently broken\n", argv [index]);
         debug_energy = filter_energies = true;
         output_format = ENERGY;
      }

      else
      {
         sprintf(comment, "unknown option `%s'", argv [index]);
         usage(comment);
         exit(102);
      }

      index++;
   }

   if (save_last_knot)
      coord_last = (ivector *) calloc(MAXSIZE, sizeof (ivector));

   //  if (filter_energies)
   //    energy_blurt_values ();

   // if no more arguments, then read from stdin instead of opening file

   if (read_from_stdin && index != argc)
   {
      sprintf(comment, "unexpected extra argument `%s'", argv [index]);
      usage(comment);
   }

   if (!read_from_stdin && index == argc)
      usage("");

   if (read_from_stdin)
   {
      fpin = stdin;
      get_a_knot = get_a_knot_binary;
      input_format = BINARY;
   }
   else
      fpin = open_file(argv [index++], "rb");

   // if no more arguments, then write to stdout instead of opening file


   if (index == argc)
   {
      fpout = stdout;
   }
   else if (index == argc - 1)
   {
      output_filename = argv [index++];
      if (numb_sequences == 0)
         fpout = open_file(output_filename, "wb");
   }
   else
   {
      sprintf(comment, "unexpected extra argument `%s'", argv [index + 1]);
      usage(comment);
   }

   // if not reading from stdin, see what kind of input we are dealing with
   if (fpin != stdin)
   {
      bool rewind_file = true;
      char magic [6];
      fread(magic, 1, 4, fpin);

      // if magic string is "CUBE" then we have binary input,
      // otherwise it is plain text
      if (!strncmp(magic, "CUBE", 4))
      {
         get_a_knot = get_a_knot_binary;
         input_format = BINARY;
      }
      else if (!strncmp(magic, "UofS", 4))
      {
         get_a_knot = get_a_knot_UofS;
         input_format = UofS;
         rewind_file = false;
         rewind(fpin);
         fgets(comment, 102, fpin); // get rid of first line
      }
      else if (magic [0] == 'N' ||
              magic [0] == 'E' ||
              magic [0] == 'W' ||
              magic [0] == 'S' ||
              magic [0] == 'U' ||
              magic [0] == 'D')
      {
         get_a_knot = get_a_knot_newsud;
         input_format = NEWSUD;
      }
      else
      {
         get_a_knot = get_a_knot_text;
         input_format = PLAIN_TEXT;
      }

      if (rewind_file)
         rewind(fpin); // go to start of input file
   }

   // if the user has not yet set the output format,
   // set the output format to be the opposite of the input format
   // UNLESS output file ends in ".b", ".bin" or ".txt"

   if (output_format == FORMAT_NOT_YET_SET)
   {
      if (ends_in(".txt", output_filename))
         output_format = PLAIN_TEXT;
      else if (ends_in(".b", output_filename) || ends_in(".bin", output_filename))
         output_format = BINARY;
      else if (input_format == PLAIN_TEXT || input_format == NEWSUD)
         output_format = BINARY;
      else
         output_format = PLAIN_TEXT;
   }

   switch (output_format)
   {
      case PLAIN_TEXT:
         put_a_knot = put_a_knot_text;
         break;
      case BINARY:
         put_a_knot = put_a_knot_binary;
         break;
      case MULTIPLE_BINARY:
         put_a_knot = put_a_knot_multiple_binary;
         break;
      case LENGTH:
         put_a_knot = put_a_knot_length;
         break;
      case ACN:
         put_a_knot = put_a_knot_ACN;
         break;
      case WRITHE:
         put_a_knot = put_a_knot_writhe;
         break;
      case ROG:
         put_a_knot = put_a_knot_rog;
         break;
         //  case ENERGY:
         //    put_a_knot = put_a_knot_energy;
         //    break;
      case BBOX:
         put_a_knot = put_a_knot_bbox;
         break;
      case PPM:
         put_a_knot = put_a_knot_PPM;
         break;
      case UofS:
         fprintf(fpout, "UofS\n"); // fall thru
      case NEWSUD:
         put_a_knot = put_a_knot_NEWSUD_or_UofS;
   }

   /*
   if (output_format == ENERGY && (flags & PROVIDE_INFO))
   fprintf (stderr, "filtering knots according to energy model %d with parameter %f\n",
   energy_model, energy_parameter);
    */

   while (get_a_knot(fpin))
      put_a_knot(fpout, flags);

   if (output_format == UofS)
      fprintf(fpout, "-999\n");

   if (fpin != stdin)
      fclose(fpin);
   if (fpout != stdout)
      fclose(fpout);

   if (numb_sequences > 0)
   {
      for (int i = 0; i < numb_sequences; i++)
      {
         fprintf(stderr, "%d ", mult_fp [i].kount);
         fclose(mult_fp [i].fp);
      }
      fprintf(stderr, "\n");
   }

   if (flags & PROVIDE_INFO)
   {
      /*
      if (output_format == ENERGY)
      fprintf (stderr, "%d knots skipped because outside energy range\n", energy_number_skipped);
       */

      if (min_length != 0 || max_length != HUGE)
         fprintf(stderr, "binning knots in length range %d <= length <= %d\n", min_length, max_length);

      fprintf(stderr, "%d knots converted from ", num_converted);
      if (input_format == PLAIN_TEXT)
         fprintf(stderr, "plain text to ");
      else if (input_format == NEWSUD)
         fprintf(stderr, "NEWSUD to ");
      else
         fprintf(stderr, "binary to ");

      switch (output_format)
      {
         case PLAIN_TEXT:
            fprintf(stderr, "plain text\n");
            break;
         case BINARY:
            fprintf(stderr, "binary\n");
            break;
         case LENGTH:
            fprintf(stderr, "lengths only\n");
            break;
         case WRITHE:
            fprintf(stderr, "writhes only\n");
            break;
         case ENERGY:
            fprintf(stderr, "energies only\n");
            break;
         case BBOX:
            fprintf(stderr, "bounding box information\n");
      }


      fprintf(stderr, "global min: (%d, %d, %d),  global max: (%d, %d, %d)\n",
              global_min [0], global_min [1], global_min [2],
              global_max [0], global_max [1], global_max [2]);
   }

   if (output_format == BBOX)
   {
      printf("volumes from %d to %d\n", bbox_min_volume, bbox_max_volume);
      int minr = HUGE, maxr = 0;
      for (int c = 0; c < 3; c++)
      {
         printf("%c-range from %d to %d\n", (char) ((int) 'x' + c), bbox_min_range [c], bbox_max_range [c]);
         minr = MIN(minr, bbox_min_range [c]);
         maxr = MAX(maxr, bbox_max_range [c]);
      }
      printf("xyz-range from %d to %d\n", minr, maxr);
   }

   if (num_tested)
      fprintf(stderr, "%d of %d knots tested were rejected (%.1f%%)\n",
           num_rejected, num_tested, 100.0 * (double) num_rejected / (double) num_tested);

   return 0;
}


