#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>

bool read_unsigned_short (unsigned short *, FILE *);
void write_unsigned_short (unsigned short *, FILE *);
bool read_int (int *, FILE *);
void write_int (int *, FILE *);

typedef struct {
  FILE *fp;
  unsigned char ch;
  int index;
} BitFile;

typedef BitFile *BitFilePtr;

int read_bit (BitFilePtr bfp);
void write_bit (int value, BitFilePtr bfp);


void panic_exit (char *s);
void panic_exit (char *s, int d);
void panic_exit (char *s, char d);
void panic_exit (char *s, char *d);
void panic_exit (char *s, char *t, int w);
bool alldigits (char *s);

FILE *open_file (char *filename, char *access);

void allocate_string (char **ps, char *what);

extern int knotnumber;

#endif

