#ifndef RANDOM_H
#define RANDOM_H

typedef double real;

extern void set_sRand_seed_to_clocktime ();
extern int rand_integer (int, int);
extern real rand_real (real low, real high);
extern real rand_uniform ();
extern void sRandSimple (int seed);

#endif
