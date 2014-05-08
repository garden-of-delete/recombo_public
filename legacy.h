#ifndef LEGACY_H
#define LEGACY_H

#include <cstdio>

// From sequence convert project types.h

typedef double real;    


typedef float  freal;   /* Always typedefed to `float'. */
typedef double dreal;   /* Always typedefed to `double'. */

typedef real vector [3];
typedef real point [3];
typedef real vector2 [2];
typedef real point2 [2];

typedef struct {
  vector v;
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

// From sequence convert project vector.h

#define PI           3.14159265358979323846
#define THREEPIBY2   4.7123889803846898576
#define TWOPI        6.28318530717958647693
#define PIBY2        1.57079632679489661923
#define PIBY4        0.78539816339744830961
#define DEG_TO_RAD   0.0174532925199432957692  /*  (PI / 180)  */
#define RAD_TO_DEG   57.2957795130823208768    /*  (180 / PI)  */
#define GOLDEN_RATIO  1.618033988749894848205

#define normal_vector(a, b)  new_normal_vector(a, b, __LINE__, __FILE__)

extern void get_sphere (vector centre, real &radius, vector *p);
extern int is_orthonormal_triad (vector u, vector v, vector w);
extern void output_triad (vector u, vector v, vector w);
extern real length_segment (vector, vector);
extern void horizon_point (vector hp, vector cent, vector eye, vector up, real radius);
extern void negate_vector (vector, vector);
extern void conv_to_xyz (vector, real, real, real);
extern void abs_vector (vector, vector);
extern void max_vector (vector, vector, vector);
extern void min_vector (vector, vector, vector);
extern void copy_vector (vector, vector);
extern void set_vector (vector, real, real, real);
extern void unit_vector (real, real, vector);
extern int equal_vector (vector, vector);
extern int equal_vector_within_epsilon (vector, vector);
extern void random_unit_vector (vector);
//extern void random_unit_vector2 (vector);
extern void random_vector (vector, real);
extern void zero_vector (vector);
extern void cross (vector, vector, vector);
extern void perp_vector (vector, vector);
extern void add_vector (vector, vector, vector);
extern void sub_vector (vector, vector, vector);
extern void mult_vector (vector, vector, real);
extern void midpoint (vector, vector, vector);
extern real dot (vector, vector);
extern real magnitude (vector);
extern void normalize_vector (vector, vector);
extern void interp_vector (vector, vector, vector, real);
extern real interp_real (real, real, real);
extern void make_orthonormal_triad (vector, vector, vector);
extern void random_orthonormal_triad (vector, vector, vector);
extern void shift_vector (vector, int);
extern void change_coord_system_inv (vector, vector, 
				     vector, vector, vector);
extern void change_coord_system (vector, vector, 
				 vector, vector, vector);

extern int is_a_zero_vector (vector);
extern int is_a_unit_vector (vector);

extern void new_normal_vector (vector, vector, int, char *);
extern void random_unit_vector_cone (vector v, real angle);
extern void test_random_unit_vector (int, int);
extern void print_vector_a (FILE *, vector);
extern void print_vector (vector);
extern void put_vector (vector);
extern void compute_normal (vector, vector, vector, vector);
extern real diff_distance_SQ (vector, vector);
extern void perp_vector_wqed (vector, vector);

extern void conv_cylindrical_to_xyz (vector, real, real, real);
extern void conv_to_cylindrical (vector, real *, real *, real *);
extern void conv_to_polar (vector, real *, real *, real *);

extern void copy_from_fvector (vector v, fvector w)    ;
extern void copy_to_fvector (fvector v, vector w);
extern void set_fvector (fvector, real, real, real);

extern vector iHat, jHat, kHat, xHat, yHat, zHat, ZERO_vect;
#define ORIGIN ZERO_vect

extern void rotate_vector_2D (vector, vector, real);
extern void set_ivector2 (ivector2, int, int);
extern void set_ivector4 (ivector4, int, int, int, int);
extern void copy_ivector4 (ivector4, ivector4);
extern void transform_point (vector, Matrix, vector);

extern int is_right_of (vector, vector, vector);
extern int inside_rect (vector, vector, real, vector);
extern void unit_vector_in_fan (vector, vector, point, real);

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
