#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#pragma warning (disable: 4996)
#pragma warning (disable: 4267)
#pragma warning (disable: 4018)
#endif

#include "id.h"

#define DEFAULT_maxlen 30000

void same_polys ();

void usage () {
  fprintf (stderr, "idknot version %s %s\n\n", __DATE__, __TIME__);
  fprintf (stderr, "Usage: idknot [options] POLYFILE > DISTRFILE\n");
  fprintf (stderr, "       idknot [options] POLYFILE DISTRFILE\n");
  fprintf (stderr, "       idknot -- [options] < POLYFILE > DISTRFILE\n");
  fprintf (stderr, "       some-poly-producing-programme | idknot -- [options] > DISTRFILE\n");
  fprintf (stderr, "       where POLYFILE is a plain text file containing a list of HOMFLY polynomials\n");
  fprintf (stderr, "       and DISTRFILE is the knot distribution\n");
  fprintf (stderr, "       and options are:\n");
  fprintf (stderr, "       -k       print knot type of each knot, line by line\n");
  fprintf (stderr, "       -ks      same as option -k above, but also print summary\n");
  fprintf (stderr, "       -u       include number of unknown knots in the summary\n");
  fprintf (stderr, "       -p       only id/count prime knots (consider composite knots as unknowns)\n");
  fprintf (stderr, "       -nz      only print non-zero counts in the summary\n");
  fprintf (stderr, "       -P       print percentages\n");
  fprintf (stderr, "       --       read from standard input instead of trying to open a data file\n");
  fprintf (stderr, "       -x KNOT  exclude (ignore) knot type KNOT\n");
  fprintf (stderr, "       -L       list all known knot types\n");
  fprintf (stderr, "       -sp      for each known knot type, print out all other knot types with same polynomial\n");
  fprintf (stderr, "\nFor full documentation and examples see http://ewok.sfsu.edu/wiki/index.php/Idknot\n");
  exit (1022);
}

void list_all_knots () {
  for (int i = 0; i < NUM_KNOTS; i++) {
    printf ("%s\n", knot [i]);
  }
  exit (111);
}

int main (int argc, char **argv) {
  int maxlen = DEFAULT_maxlen;
  bool printeach = false;
  bool printsummary = true;
  bool printunknown = false;
  bool non_zero_only = false;
  bool read_from_stdin = false;
  bool print_percent = false;

  int nknots = NUM_KNOTS;
  char the_knot [22];
  int num_unknown = 0;
  char *exclude = (char *) NULL;

  int index = 1;

  while (index < argc && argv [index][0] == '-') {
    if (!strcmp (argv [index], "-k")) {
      printsummary = false;
      printeach = true;
    }
    else if (!strcmp (argv [index], "-ks")) {
      printsummary = true;
      printeach = true;
    }
    else if (!strcmp (argv [index], "-u")) 
      printunknown = true;
    else if (!strcmp (argv [index], "-P")) 
      print_percent = true;
    else if (!strcmp (argv [index], "--")) 
      read_from_stdin = true;
    else if (!strcmp (argv [index], "-nz")) 
      non_zero_only = true;
    else if (!strcmp (argv [index], "-p")) 
      nknots = NUM_PRIME_KNOTS;
    else if (!strcmp (argv [index], "-x")) {
      if (index + 1 < argc) 
	exclude = argv [++index];
      else 
	usage ();
    }
    else if (!strcmp (argv [index], "-L"))
      list_all_knots ();
    else if (!strcmp (argv [index], "-sp"))
      same_polys ();
    else {
      fprintf (stderr, "unknown option `%s'\n", argv [index]);
      exit (1022);
    }

    index++;
  }

  FILE *fpout;
  FILE *fp;

  if (read_from_stdin && index != argc || !read_from_stdin && index == argc) usage ();

  if (read_from_stdin) 
    fp = stdin;
  else {
    fp = fopen (argv [index], "r");
    if (!fp) {
      fprintf (stderr, "Failed to open file `%s' for read access\n", argv [index]);
      exit (108);
    }
    index++;
  }


  if (index == argc) {
    fpout = stdout;
  }
  else if (index == argc - 1) {
    fpout = fopen (argv [index], "w");
    if (!fpout) {
      fprintf (stderr, "Failed to open file `%s' for write access\n", argv [index]);
      exit (108);
    }
  }
  else {
    usage ();
    exit (102);
  }


  char *inp = (char *) calloc (maxlen + 2, sizeof (char));
  int *kount = (int *) calloc (NUM_KNOTS, sizeof (int));
  int total = 0;

  int excluded_ID = NOTHING;   // note that NOTHING is different from UNKNOWN
  if (exclude) {
    for (int i = 0; i < NUM_KNOTS; i++) {
      if (!strcmp (exclude, knot [i])) {
	excluded_ID = i;
	break;
      }
    }
  }

  while (fgets (inp, maxlen, fp)) {
    // strip off the new line character at end of string
    inp [strlen (inp) - 1] = (char) 0;

    int knotID = id_knot (the_knot, inp, nknots);

    if (knotID == excluded_ID) continue;

    if (knotID == UNKNOWN) 
      ++num_unknown;
    else
      kount [knotID]++;
    
    ++total;
    if (printeach)
      fprintf (fpout, "%s\n", the_knot);
  }

  if (printsummary) {
    for (int i = 0; i < nknots; i++) {
      if (!non_zero_only || kount [i] != 0) {
	fprintf (fpout, "%s %d", knot [i], kount [i]);
	if (print_percent) 
	  fprintf (fpout, " %6.2f%%\n", (double) kount [i] * 100.0 / (double) total);
	else
	  fprintf (fpout, "\n");
      }
    }
    if (printunknown && (!non_zero_only || num_unknown != 0)) {
      fprintf (fpout, "unknown %d", num_unknown);
      if (print_percent) 
	fprintf (fpout, " %6.2f%%\n", (double) num_unknown * 100.0 / (double) total);
      else
	fprintf (fpout, "\n");
    }
    fprintf (fpout, "total %d", total);
    if (print_percent) 
      fprintf (fpout, " %6.2f%%\n", (double) total * 100.0 / (double) total);
    else
      fprintf (fpout, "\n");
  }

  fclose (fp);
  if (fpout != stdout)
    fclose (fpout);

  return 0;
}
