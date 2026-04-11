#ifndef LEGACY_H
#define LEGACY_H

#include <cstdio>

// From sequence convert project types.h

typedef double real;

typedef int ivector [3];

typedef real vector3 [3];

// From sequence convert project vector3.h

#define PI           3.14159265358979323846
#define THREEPIBY2   4.7123889803846898576
#define TWOPI        6.28318530717958647693
#define PIBY2        1.57079632679489661923
#define PIBY4        0.78539816339744830961
#define DEG_TO_RAD   0.0174532925199432957692  /*  (PI / 180)  */
#define RAD_TO_DEG   57.2957795130823208768    /*  (180 / PI)  */
#define GOLDEN_RATIO  1.618033988749894848205

extern void cross (vector3, vector3, vector3);
extern void sub_vector (vector3, vector3, vector3);
extern real dot (vector3, vector3);
extern real magnitude (vector3);
extern void normalize_vector (vector3, vector3);
extern void set_vector (vector3, real, real, real);

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

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))

#endif
