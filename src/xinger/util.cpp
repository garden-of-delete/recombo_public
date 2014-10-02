#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef _WIN32   // disable annoying warnings
#pragma warning (disable: 4996)
#pragma warning (disable: 4101)
#pragma warning (disable: 4267)
#endif

#include "util.h"

int knotnumber = 1;

bool read_unsigned_short (unsigned short *v, FILE *fp) {
  unsigned char uca, ucb;
  if (fread (&ucb, 1, 1, fp) != 1) return false;
  if (fread (&uca, 1, 1, fp) != 1) return false;
  *v = (ucb << 8) | uca;
  return true;
}

void write_unsigned_short (unsigned short *val, FILE *fp) {
  unsigned char a, b;
  a = (unsigned char) (0xff & (*val) >> 8);
  b = (unsigned char) (0xff & (*val));

  fwrite (&a, 1, 1, fp);
  fwrite (&b, 1, 1, fp);
}

void write_int (int *val, FILE *fp) {
  unsigned char a, b, c, d;

  a = (unsigned char) (0xff & (*val) >> 24);
  b = (unsigned char) (0xff & (*val) >> 16);
  c = (unsigned char) (0xff & (*val) >> 8);
  d = (unsigned char) (0xff & (*val));

  fwrite (&a, 1, 1, fp);
  fwrite (&b, 1, 1, fp);
  fwrite (&c, 1, 1, fp);
  fwrite (&d, 1, 1, fp);
}

bool read_int (int *v, FILE *fp) {
  unsigned char uc;
  *v = 0;
  for (int shift = 24; shift >= 0; shift -= 8) {
    if (fread (&uc, 1, 1, fp) != 1) return false;
    *v |= uc << shift;
  }
  return true;
}

void panic_exit (char *s) {
  fprintf (stderr, "\n *** %s (was reading knot #%d)\n\n", s, knotnumber);
  fflush (stderr);
  exit (102);
}

void panic_exit (char *s, int d) {
  char *t = (char *) calloc (strlen (s) + 108, sizeof (char));
  sprintf (t, s, d);
  panic_exit (t);
}

void panic_exit (char *r, char *s, int d) {
  char *t = (char *) calloc (strlen (s) + strlen (r) + 108, sizeof (char));
  sprintf (t, r, s, d);
  panic_exit (t);
}

void panic_exit (char *s, char d) {
  char *t = (char *) calloc (strlen (s) + 108, sizeof (char));
  sprintf (t, s, d);
  panic_exit (t);
}

void panic_exit (char *s, char *d) {
  char *t = (char *) calloc (strlen (s) + strlen (d) + 108, sizeof (char));
  sprintf (t, s, d);
  panic_exit (t);
}

FILE *open_file (char *filename, char *access) {
  FILE *fp = fopen (filename, access);
  if (!fp) {
    fprintf (stderr, " *** Can't open file `%s' for %s access!\n", filename, access);
    fflush (stderr);
    exit (1);
  }

  return fp;
}


void allocate_string (char **ps, char *what) {
  if (*ps)
    free (*ps);
  *ps = (char *) NULL;
  if (what) {
    *ps = (char *) calloc (strlen (what) + 2, sizeof (char));
    strcpy (*ps, what);
  }
}


void write_bit (int b, BitFilePtr bfp) {
  bfp->ch |= (b & 0x1) << bfp->index;
  if (++bfp->index == 8) {
    //fprintf (stderr, "wrote %d\n", bfp->ch);
    fwrite (&bfp->ch, 1, 1, bfp->fp);
    bfp->ch = 0;
    bfp->index = 0;
  }
}

int read_bit (BitFilePtr bfp) {
  if (bfp->index % 8 == 0) {
    if (fread (&bfp->ch, 1, 1, bfp->fp) != 1) panic_exit ("end of file before end of data");
    bfp->index = 0;
  }
  
  int value = bfp->ch >> bfp->index;
  bfp->index++;

  return 0x1 & value;
}

bool alldigits (char *s) {
  while (*s) {
    if (!isdigit ((int) *s)) return false;
    s++;
  }
  return true;
}
