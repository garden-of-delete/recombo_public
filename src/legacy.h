#ifndef LEGACY_H
#define LEGACY_H

#include <cstdio>

// From sequence convert project types.h

typedef double real;    


typedef float  freal;   /* Always typedefed to `float'. */
typedef double dreal;   /* Always typedefed to `double'. */

typedef real vector3 [3];
typedef real point [3];
typedef real vector2 [2];
typedef real point2 [2];

typedef struct {
  vector3 v;
  real w;
} vector4;

typedef unsigned short usvector [3];
typedef  short svector [3];
typedef float  fvector [3];
typedef float  fvector2 [2];
typedef double dvector [3];
typedef real Mat3x3 [9];
typedef signed char scvector [3];

typedef int ivector4 [4];
typedef int ivector [3];
typedef int ivector2 [2];

typedef struct {
  int numer;
  int denom;
} rational;

typedef real Matrix [16];

typedef char *TextString;
typedef TextString *TextStringPtr;

// From sequence convert project vector3.h

#define PI           3.14159265358979323846
#define THREEPIBY2   4.7123889803846898576
#define TWOPI        6.28318530717958647693
#define PIBY2        1.57079632679489661923
#define PIBY4        0.78539816339744830961
#define DEG_TO_RAD   0.0174532925199432957692  /*  (PI / 180)  */
#define RAD_TO_DEG   57.2957795130823208768    /*  (180 / PI)  */
#define GOLDEN_RATIO  1.618033988749894848205

#define normal_vector(a, b)  new_normal_vector(a, b, __LINE__, __FILE__)

extern void get_sphere (vector3 centre, real &radius, vector3 *p);
extern int is_orthonormal_triad (vector3 u, vector3 v, vector3 w);
extern void output_triad (vector3 u, vector3 v, vector3 w);
extern real length_segment (vector3, vector3);
extern void horizon_point (vector3 hp, vector3 cent, vector3 eye, vector3 up, real radius);
extern void negate_vector (vector3, vector3);
extern void conv_to_xyz (vector3, real, real, real);
extern void abs_vector (vector3, vector3);
extern void max_vector (vector3, vector3, vector3);
extern void min_vector (vector3, vector3, vector3);
extern void copy_vector (vector3, vector3);
extern void set_vector (vector3, real, real, real);
extern void unit_vector (real, real, vector3);
extern int equal_vector (vector3, vector3);
extern int equal_vector_within_epsilon (vector3, vector3);
extern void random_unit_vector (vector3);
//extern void random_unit_vector2 (vector3);
extern void random_vector (vector3, real);
extern void zero_vector (vector3);
extern void cross (vector3, vector3, vector3);
extern void perp_vector (vector3, vector3);
extern void add_vector (vector3, vector3, vector3);
extern void sub_vector (vector3, vector3, vector3);
extern void mult_vector (vector3, vector3, real);
extern void midpoint (vector3, vector3, vector3);
extern real dot (vector3, vector3);
extern real magnitude (vector3);
extern void normalize_vector (vector3, vector3);
extern void interp_vector (vector3, vector3, vector3, real);
extern real interp_real (real, real, real);
extern void make_orthonormal_triad (vector3, vector3, vector3);
extern void random_orthonormal_triad (vector3, vector3, vector3);
extern void shift_vector (vector3, int);
extern void change_coord_system_inv (vector3, vector3, 
				     vector3, vector3, vector3);
extern void change_coord_system (vector3, vector3, 
				 vector3, vector3, vector3);

extern int is_a_zero_vector (vector3);
extern int is_a_unit_vector (vector3);

extern void new_normal_vector (vector3, vector3, int, char *);
extern void random_unit_vector_cone (vector3 v, real angle);
extern void test_random_unit_vector (int, int);
extern void print_vector_a (FILE *, vector3);
extern void print_vector (vector3);
extern void put_vector (vector3);
extern void compute_normal (vector3, vector3, vector3, vector3);
extern real diff_distance_SQ (vector3, vector3);
extern void perp_vector_wqed (vector3, vector3);

extern void conv_cylindrical_to_xyz (vector3, real, real, real);
extern void conv_to_cylindrical (vector3, real *, real *, real *);
extern void conv_to_polar (vector3, real *, real *, real *);

extern void copy_from_fvector (vector3 v, fvector w)    ;
extern void copy_to_fvector (fvector v, vector3 w);
extern void set_fvector (fvector, real, real, real);

extern vector3 iHat, jHat, kHat, xHat, yHat, zHat, ZERO_vect;
#define ORIGIN ZERO_vect

extern void rotate_vector_2D (vector3, vector3, real);
extern void set_ivector2 (ivector2, int, int);
extern void set_ivector4 (ivector4, int, int, int, int);
extern void copy_ivector4 (ivector4, ivector4);
extern void transform_point (vector3, Matrix, vector3);

extern int is_right_of (vector3, vector3, vector3);
extern int inside_rect (vector3, vector3, real, vector3);
extern void unit_vector_in_fan (vector3, vector3, point, real);

bool equal_ivector (ivector a, ivector b);
void sub_ivector (ivector v, ivector w, ivector z);
void add_ivector (ivector v, ivector w, ivector z);
void copy_ivector (ivector v, ivector w);
void negate_ivector (ivector v, ivector w);

void mult_ivector (ivector v, ivector w, double con);
void mult_ivector (ivector v, ivector w, int con);
void midpoint (ivector mid, ivector a, ivector b);
void max_ivector (ivector v, ivector w, ivector z);
void min_ivector (ivector v, ivector w, ivector z);
void set_ivector (ivector a, int x, int y, int z);


extern void mat3X3mult (Mat3x3, Mat3x3, Mat3x3);
extern void mat3X3inv (Mat3x3);
extern real mat3x3det (Mat3x3 mat);

// From sequence convert project util.h

typedef struct {
  FILE *fp;
  unsigned char ch;
  int index;
} BitFile;

typedef BitFile *BitFilePtr;

#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

extern void read_unsigned_short (unsigned short *, FILE *);
extern void read_int (int *, FILE *);
extern void write_int (int *, FILE *);
extern void read_float (float *, FILE *);
extern void read_double (double *, FILE *);
extern void read_fvector (fvector, FILE *);
extern void write_ivector (ivector v, FILE *fp);
extern void read_ivector (ivector v, FILE *fp);
extern void copy_ivector (ivector a, ivector b);

extern bool output_CUBE (int num_vertices, FILE *fp, ivector *coord);
extern bool input_CUBE (int &num_vertices, FILE *fp, ivector *coord);

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))


// From sequence convert project writhe.h

void bbox (int &volume, ivector range, int nvert, ivector *vert);
double writhe (double &ACN, int nvert, ivector *vert, double jitter);
double radius_of_gyration (int nvert, ivector *vert);

#endif
