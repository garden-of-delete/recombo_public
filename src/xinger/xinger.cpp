#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32   // disable annoying warnings
#pragma warning (disable: 4996)
#pragma warning (disable: 4101)
#endif

#include "xinger.h"
#include "random.h"
#include "util.h"

bool support_links = false;
int number_of_components;
XingPtr comp_xp [10];
int     comp_strand [10];

FILE *fpkps = (FILE *) NULL;
FILE *fp_ray = (FILE *) NULL;

#define RI  1
#define RII 2
#define RRC 3
int which_reid_move = 0;
void removeXings (DiagramPtr diag, bool R2);

bool show_iso_clasps = true;
#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

int nmn (int i, int n) {
  // returns the next 
  int j = abs (i) - 1;  // convert to zero relative
  return (j + 1) % n + 1;
}

int pmn (int i, int n) {
  // returns the prev 
  int j = abs (i) - 1;  // convert to zero relative
  return (j + n - 1) % n + 1;
}
bool superv = false;

int abovebelow (XingPtr xp, int strand) {
  // returns 1 is crossing is an overcrossing
  // and -1 if an undercrossing

  if (xp->sign == 1) {
    if (strand == 0) 
      return 1;
    else 
      return -1;
  }
  if (strand == 1) 
    return 1;
  else 
    return -1;
}


// these routines may be used in the fashion of
// fore (xp, strand, xp, strand);
// without causing any problems

void fore (XingPtr xpi, int strand, XingPtr &fxp, int &fstrand) {
  XingPtr xp = xpi;
  fxp = xp->xtop [strand];
  fstrand = xp->stop [strand];
}

void back (XingPtr xpi, int strand, XingPtr &bxp, int &bstrand) {
  XingPtr xp = xpi;
  bxp = xp->xbot [strand];
  bstrand = xp->sbot [strand];
}

void right (XingPtr xpi, int strand, XingPtr &rxp, int &rstrand) {
  XingPtr xp = xpi;
  if (strand == 0) {
    rxp = xp->xbot [1];
    rstrand = xp->sbot [1];
  }
  else {
    rxp = xp->xtop [0];
    rstrand = xp->stop [0];
  }
}

void left (XingPtr xpi, int strand, XingPtr &lxp, int &lstrand) {
  XingPtr xp = xpi;
  if (strand == 0) {
    lxp = xp->xtop [1];
    lstrand = xp->stop [1];
  }
  else {
    lxp = xp->xbot [0];
    lstrand = xp->sbot [0];
  }
}


DiagramPtr initDiagram (int max_xings) {
  DiagramPtr diag = (DiagramPtr) calloc (1, sizeof (Diagram));
  diag->max_xings = max_xings;
  diag->xing = (XingPtr *) calloc (diag->max_xings, sizeof (XingPtr));
  if (max_xings > 999) 
    panic_exit ("maximum number of crossings exceeded");

  diag->gauss_code = (char *) calloc (5 * max_xings + 2, sizeof (char));
  diag->cscratch = (char *) calloc (5 * max_xings + 2, sizeof (char));
  for (int i = 0; i < diag->max_xings; i++) 
    diag->xing [i] = (XingPtr) calloc (1, sizeof (Xing));

  diag->loc = (int *) calloc (2 * max_xings + 1, sizeof (int));
  diag->scratch = (int *) calloc (2 * max_xings + 1, sizeof (int));
  //fprintf (stderr, "Diagram initialized with %d max_xings\n", diag->max_xings); fflush (stderr);
  diag->p = 1.08;
  return diag;
}

DiagramPtr initDiagram (int max_xings, double p, bool zok) {
  DiagramPtr diag = initDiagram (max_xings);
  diag->p = p;
  diag->zok = zok;
  return diag;
}

void tie_ends (XingPtr xpa, int stranda, XingPtr xpb, int strandb) {
  // tie stranda of xpa (a top strand) to strandb of xpb (a bot strand)
  //printf ("tying top strand %d of Xing %d to bot strand %d of Xing %d\n", stranda, xpa->ID, strandb, xpb->ID);
  xpa->xtop [stranda] = xpb;
  xpa->stop [stranda] = strandb;
  xpb->xbot [strandb] = xpa;
  xpb->sbot [strandb] = stranda;
}

void tie_xings (XingPtr xa, int labela, XingPtr xb, int labelb) {
  // tie top of xa to bot of xb
  int stranda, strandb;

  //printf ("tying %d/%d to %d/%d\n", xa->ID, labela, xb->ID, labelb);
  if (!xa) panic_exit ("NULL A");
  if (!xb) panic_exit ("NULL B");

  if (!(xb == xa->xtop [0] || xb == xa->xtop [1]))
    panic_exit ("tie_xings() failure A");

  if (!(xa == xb->xbot [0] || xa == xb->xbot [1]))
    panic_exit ("tie_xings() failure B");

  if (labela == xa->lab [0])
    stranda = 0;
  else if (labela == xa->lab [1])
    stranda = 1;
  else 
    panic_exit ("tie_xings() failure label A");

  if (labelb == xb->lab [0])
    strandb = 0;
  else if (labelb == xb->lab [1])
    strandb = 1;
  else 
    panic_exit ("tie_xings() failure label B");

  xa->stop [stranda] = strandb;
  xb->sbot [strandb] = stranda;
}

#define SWAP(A,B)  {tmp = A; A = B; B = tmp;}

void makeDiagramFromGauss (DiagramPtr diag, int *gauss, int *sign, Coord *coords, int nxings) {
  makeDiagramFromGauss (diag, gauss, sign, nxings);
  int index = 0;
  int n2 = 2 * nxings;
  for (int i = 1; i <= n2; i += 2, index++) {
    XingPtr xp = diag->xing [index];    // convenience pointer
    xp->x = coords [i].x;
    xp->y = coords [i].y;
    xp->z = coords [i].z;
    xp->angle = coords [i].angle;
    //fprintf (stderr, "Xing #%d located at (%f, %f) and angle %f\n", xp->ID, xp->x, xp->y, xp->angle);
  }
}

void makeDiagramFromGauss (DiagramPtr diag, int *gauss, int *sign, int nxings) {
  diag->nxings = nxings;
  diag->first_xing = diag->xing [0];
  int index = 0;
  int n2 = 2 * nxings;
  for (int i = 1; i <= n2; i += 2, index++) {
    XingPtr xp = diag->xing [index];    // convenience pointer

    xp->ID = index;
    xp->location = index;
    xp->lab [0] = xp->lab [1] = xp->stop [0] = xp->stop [1] = xp->sbot [0] = xp->sbot [1] = xp->start = 108;

    xp->sign = sign [i];

    int strand0 = i;
    int strand1 = abs (gauss [i]);
    int tmp;

    if (gauss [i] < 0) SWAP (strand0, strand1);
    if (xp->sign == -1) SWAP (strand0, strand1);

    xp->lab [0] = strand0;
    xp->lab [1] = strand1;

    if (strand0 < strand1)
      xp->start = 0;
    else
      xp->start = 1;
    
    diag->loc [i] = diag->loc [abs (gauss [i])] = index;
  }
  if (index != diag->nxings) panic_exit ("complete failure!");

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];    // convenience pointer
    

    xp->xtop [0] = diag->xing [diag->loc [nmn (xp->lab [0], n2)]];
    xp->xtop [1] = diag->xing [diag->loc [nmn (xp->lab [1], n2)]];
    xp->xbot [0] = diag->xing [diag->loc [pmn (xp->lab [0], n2)]];
    xp->xbot [1] = diag->xing [diag->loc [pmn (xp->lab [1], n2)]];
  }
  for (int i = 1; i <= n2; i++) {
    int inext = nmn (i, n2);
    XingPtr xpa = diag->xing [diag->loc [i]];
    XingPtr xpb = diag->xing [diag->loc [inext]];
    tie_xings (xpa, i, xpb, inext);
    gauss [i - 1] = gauss [i] = sign [i - 1] = sign [i] = 0;   // clear things for next knot 
  }
}

void mark_all (DiagramPtr diag, int value) {
  for (int i = 0; i < diag->nxings; i++) 
    diag->xing [i]->mark = value;
}

#define IMPOSSIBLE_LABEL 16008

bool labelDiagramLink (DiagramPtr diag) {
  XingPtr xp, xpn;

  for (int i = 0; i < diag->nxings; i++) {  // clobber the old labels
    xp = diag->xing [i];
    xp->lab [0] = xp->lab [1] = IMPOSSIBLE_LABEL;
    xp->gauss = 0;
    xp->mark = 0;
    xp->test = 0;
    xp->last_strand = 0;
  }
  
  XingPtr xpstart;

  xp = xpstart = comp_xp [0] = diag->first_xing;
  int label = 1;
  int strand = xp->start;
  comp_strand [0] = strand;
  int gauss = 1;
 
  number_of_components = 1;

  while (label <= 2 * diag->nxings) {
    xp->test++;
    xp->lab [strand] = label;
    xpn = xp->xtop [strand];
    xp->last_strand = strand;
    strand = xp->stop [strand];
    if (xp->gauss == 0)
      xp->gauss = gauss++;   // only assign and increment first time Xing is encountered

    if (xpn == comp_xp [number_of_components - 1] && strand == comp_strand [number_of_components - 1]) {
      for (int isis = 0; isis < diag->nxings; isis++) {
	if (diag->xing [isis]->test < 2) {
	  xp = xpstart = comp_xp [number_of_components] = diag->xing [isis];
	  strand = comp_strand [number_of_components] = 1 - xp->last_strand;
	  number_of_components++;
	  break;
	}
      }
    }
    else
      xp = xpn;
 
    label++;
  }

  // check labels
  for (int i = 0; i < diag->nxings; i++) {
    if (diag->xing [i]->lab [0] == IMPOSSIBLE_LABEL) 
      panic_exit ("Xing %d has impossible lab [0]", i);
    if (diag->xing [i]->lab [1] == IMPOSSIBLE_LABEL) 
      panic_exit ("Xing %d has impossible lab [1]", i);
    if (diag->xing [i]->gauss == 0) 
      panic_exit ("Xing %d has zero gauss label", i);
    if (diag->xing [i]->test != 2) 
      panic_exit ("Xing %d has wrong visit count", i);
    //printf ("%d:  %d %d\n", i, diag->xing [i]->lab [0],  diag->xing [i]->lab [1]);
  }
  return true;
}

bool labelDiagram (DiagramPtr diag) {
  if (support_links) {
    return labelDiagramLink (diag);
  }
  XingPtr xp, xpn;

  for (int i = 0; i < diag->nxings; i++) {  // clobber the old labels
    xp = diag->xing [i];
    xp->lab [0] = xp->lab [1] = IMPOSSIBLE_LABEL;
    xp->gauss = 0;
    xp->test = 0;
  }

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;
  int gauss = 1;

  while (label <= 2 * diag->nxings) {
    xp->test++;
    xp->lab [strand] = label;
    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    if (xp->gauss == 0)
      xp->gauss = gauss++;   // only assign and increment first time Xing is encountered

    xp = xpn;
    label++;
  }

  // check labels
  for (int i = 0; i < diag->nxings; i++) {
    if (diag->xing [i]->lab [0] == IMPOSSIBLE_LABEL) 
      panic_exit ("Xing %d has impossible lab [0]", i);
    if (diag->xing [i]->lab [1] == IMPOSSIBLE_LABEL) 
      panic_exit ("Xing %d has impossible lab [1]", i);
    if (diag->xing [i]->gauss == 0) 
      panic_exit ("Xing %d has zero gauss label", i);
    if (diag->xing [i]->test != 2) 
      panic_exit ("Xing %d has wrong visit count", i);
    //printf ("%d:  %d %d\n", i, diag->xing [i]->lab [0],  diag->xing [i]->lab [1]);
  }
  return true;
}


bool labelDiagram2_not_used (DiagramPtr diag) {
  XingPtr xp, xpn;

  for (int i = 0; i < diag->nxings; i++) {  // clobber the old labels
    xp = diag->xing [i];
    xp->lab [0] = xp->lab [1] = IMPOSSIBLE_LABEL;
    xp->gauss = 0;
    xp->test = 0;
  }

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;
  int gauss = 1;

  while (label <= 2 * diag->nxings) {
    xp->test++;
    xp->lab [strand] = label;
    xpn = xp->xbot [strand];
    strand = xp->sbot [strand];
    if (xp->gauss == 0)
      xp->gauss = gauss++;   // only assign and increment first time Xing is encountered

    xp = xpn;
    label++;
  }

  // check labels
  for (int i = 0; i < diag->nxings; i++) {
    if (diag->xing [i]->lab [0] == IMPOSSIBLE_LABEL) 
      panic_exit ("2 Xing %d has impossible lab [0]", i);
    if (diag->xing [i]->lab [1] == IMPOSSIBLE_LABEL) 
      panic_exit ("2 Xing %d has impossible lab [1]", i);
    if (diag->xing [i]->gauss == 0) 
      panic_exit ("2 Xing %d has zero gauss label", i);
    if (diag->xing [i]->test != 2) 
      panic_exit ("2 Xing %d has wrong visit count", i);
  }
  return true;
}

void DTcode (FILE *fpc, FILE *fpd, FILE *fpo, DiagramPtr diag) {
  labelDiagram (diag);

  if (fpc)
    fprintf (fpc, "%d", diag->nxings);

  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;

  while (label <= 2 * diag->nxings) {
    if (label % 2 == 1) {
      if (strand == 0) {
        if (xp->sign == -1) fprintf (fpd, "-");
      }
      else {
        if (xp->sign == 1) fprintf (fpd, "-");
      }
      fprintf (fpd, "%d ", xp->lab [1 - strand]);
      if (fpo)
        fprintf (fpo, "%d ", xp->sign);
    }
    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  fprintf (fpd, "\n");
  if (fpo)
    fprintf (fpo, "\n");
}


#define BOUT_START  1
#define BOUT_XING   2
#define BOUT_FINISH 3


void DTcode_binary_output (FILE *fp, int what, int a, int b, int c) {
  static int num_bits_per_xing;
  int magic_XING = MAGIC_XING;
  unsigned char XING_version = (char) 0;
  static bool first_call = true;
  static BitFile bf;
  static unsigned short nxings;
  int max_label = 4;

  switch (what) {
  case BOUT_START:
    if (first_call) {
      first_call = false;
      write_int (&magic_XING, fp);      // magic number to indicate this is a binary XING file
      fwrite (&XING_version, 1, 1, fp);
    }
    nxings = (unsigned short) a;
    write_unsigned_short (&nxings, fp);
    bf.ch = (char) 0;
    bf.index = 0;
    bf.fp = fp;
    num_bits_per_xing = 1;
    while (2 * nxings > max_label){
      num_bits_per_xing++;
      max_label *= 2;
    }

    break;
  case BOUT_XING:
    // write sign (orientation) of Xing
    if (a == 1) 
      write_bit (1, &bf);
    else
      write_bit (0, &bf);

    // write sign of DT code
    if (b == 1)
      write_bit (1, &bf);
    else
      write_bit (0, &bf);

    if (c % 2 != 0) 
      panic_exit ("DT code label not even as expected");

    // in the binary file, a num_bits_per_xing-bit binary number is used
    // to represent the DT code label, for example, the binary numbers
    // 000 001 010 011 100 ... represent the decimal DT code labels
    // 2   4   6   8   10  ... 

    c = c / 2 - 1;
    for (int i = 0; i < num_bits_per_xing; i++) {
      write_bit (0x01 & c, &bf);
      c = c >> 1;
    }
    break;
  case BOUT_FINISH:
    int num_bits_written = nxings * (2 + num_bits_per_xing);

    // fill out the rest of the byte with 1s, if necessary
    while (num_bits_written % 8 != 0) {
      write_bit (1, &bf);
      num_bits_written++;
    }
  }
}

void DTcode (FILE *fpb, DiagramPtr diag) {  // binary output
  if (!labelDiagram (diag)) panic_exit ("diagram couldn't be labelled");

  DTcode_binary_output (fpb, BOUT_START, diag->nxings, 0, 0);

  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;

  while (label <= 2 * diag->nxings) {
    int dt_sign = 1;

    if (label % 2 == 1) {
      if (strand == 0) {
        if (xp->sign == -1) dt_sign = -1;
      }
      else {
        if (xp->sign == 1) dt_sign = -1;
      }
      DTcode_binary_output (fpb, BOUT_XING, xp->sign, dt_sign, xp->lab [1 - strand]);
    }
    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  DTcode_binary_output (fpb, BOUT_FINISH, 0, 0, 0);
}

void GaussCodeLink (FILE *fp, DiagramPtr diag) {
  if (!labelDiagramLink (diag)) panic_exit ("GaussCodeLink(): diagram couldn't be labelled");
  XingPtr xp, xpn, xpstart;
  /*
  static bool fcall = true;
  if (fcall) {
    fprintf (stderr, "GaussCodeLink()\n");
    fcall = false;
  }
  */
  xp = xpstart = diag->first_xing;
  int strand = xp->start;

  for (int comp = 0; comp < number_of_components; comp++) {
    xp = comp_xp [comp];
    strand = comp_strand [comp];

    bool looping = true;

    while (looping) {
      if (strand == 0) {
	if (xp->sign == 1) 
	  fprintf (fp, "a");
	else 
	  fprintf (fp, "b"); 
      }
      else {
	if (xp->sign == 1) 
	  fprintf (fp, "b");
	else 
	  fprintf (fp, "a");
      }
      
      fprintf (fp, "%d", xp->gauss);
      
      if (xp->sign == 1)
	fprintf (fp, "+");
      else
	fprintf (fp, "-");
      
      xpn = xp->xtop [strand];
      strand = xp->stop [strand];
      xp = xpn;
      if (xp == comp_xp [comp] && strand == comp_strand [comp]) {
	looping = false;
	if (comp < number_of_components - 1) {
	  fprintf (fp, "|");
	  xp = comp_xp [comp + 1];
	  strand = comp_strand [comp + 1];
	}
      }
    } 
  }
  fprintf (fp, "\n");
}

void GaussCode_old (FILE *fp, DiagramPtr diag) {
  if (support_links) {
    GaussCodeLink (fp, diag);
    return;
  }
  if (!labelDiagram (diag)) panic_exit ("GaussCode(): diagram couldn't be labelled");
  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;

  while (label <= 2 * diag->nxings) {
    if (strand == 0) {
      if (xp->sign == 1) 
        fprintf (fp, "a");
      else 
        fprintf (fp, "b"); 
    }
    else {
      if (xp->sign == 1) 
        fprintf (fp, "b");
      else 
        fprintf (fp, "a");
    }

    fprintf (fp, "%d", xp->gauss);

    if (xp->sign == 1)
      fprintf (fp, "+");
    else
      fprintf (fp, "-");

    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  fprintf (fp, "\n");
}

void GaussCode (FILE *fp, DiagramPtr diag) {
  if (support_links) {
    GaussCodeLink (fp, diag);
    return;
  }
  if (!labelDiagram (diag)) panic_exit ("GaussCode(): diagram couldn't be labelled");
  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;
  int index = 0;
  char buff [8];

  while (label <= 2 * diag->nxings) {
    if (strand == 0) {
      if (xp->sign == 1) 
	diag->gauss_code [index++] = 'a';
      else 
	diag->gauss_code [index++] = 'b';
    }
    else {
      if (xp->sign == 1) 
	diag->gauss_code [index++] = 'b';
      else 
	diag->gauss_code [index++] = 'a';
    }

    sprintf (buff, "%d", xp->gauss);
    int len = strlen (buff);
    for (int u = 0; u < len; u++) 
      diag->gauss_code [index++] = buff [u];
    
    if (xp->sign == 1)
      diag->gauss_code [index++] = '+';
    else
      diag->gauss_code [index++] = '-';

    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  diag->gauss_code [index] = (char) 0;
  if (fp) fprintf (fp, "%s\n", diag->gauss_code);
}

void GaussCode_recomb (FILE *fp, DiagramPtr diag) {
  if (!labelDiagram (diag)) panic_exit ("GaussCode(): diagram couldn't be labelled");
  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;

  // pick a xing at random and make it into an infinity tangle
  XingPtr inf = diag->xing [rand_integer (0, diag->nxings)];

  mark_all (diag, 0);
  
  while (label <= 2 * diag->nxings) {
    if (strand == 0) {
      if (xp->sign == 1) 
        fprintf (fp, "a");
      else 
        fprintf (fp, "b"); 
    }
    else {
      if (xp->sign == 1) 
        fprintf (fp, "b");
      else 
        fprintf (fp, "a");
    }

    fprintf (fp, "%d", xp->gauss);

    if (xp->sign == 1)
      fprintf (fp, "+");
    else
      fprintf (fp, "-");

    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  fprintf (fp, "\n");
}

bool Jenkins_old_style = true;

void JenkinsHOMFLY (FILE *fp, DiagramPtr diag) {
  if (!labelDiagram (diag)) panic_exit ("JenkinsHOMFLY(): diagram couldn't be labelled");
  XingPtr xp, xpn;

  xp = diag->first_xing;
  int label = 1;
  int strand = xp->start;

  fprintf (fp, "1\n%d\n", 2 * diag->nxings);

  while (label <= 2 * diag->nxings) {
    fprintf (fp, "%d ", xp->gauss - 1);
    if (strand == 0) {
      if (xp->sign == 1) 
        fprintf (fp, "1 ");
      else 
        fprintf (fp, "-1 "); 
    }
    else {
      if (xp->sign == 1) 
        fprintf (fp, "-1 ");
      else 
        fprintf (fp, "1 ");
    }



    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
    label++;
  } 
  if (Jenkins_old_style)
    fprintf (fp, "\n%d\n", diag->nxings);

  for (int i = 0; i < diag->nxings; i++) 
    fprintf (fp, "%d %d\n", diag->xing [i]->gauss - 1, diag->xing [i]->sign);
}

void DTcode (DiagramPtr diag) {
  DTcode ((FILE *) NULL, stdout, (FILE *) NULL, diag);
}

void makeGaussFromDiagram (int *gauss, int *sign, int nxings, DiagramPtr diag) {
  //DTcode (stdout, diag);
}

void makeGaussFromDiagram2 (int *gauss, int *sign, int nxings, DiagramPtr diag) {
  //DTcode (stdout, diag);
}

void swapXings (DiagramPtr diag, int A, int B) {
  if (A == B) return;
  XingPtr tmp;

  SWAP (diag->xing [A], diag->xing [B]);
  diag->xing [A]->location = A;
  diag->xing [B]->location = B;
}

void recombXing_notused (DiagramPtr diag) {
  if (diag->nxings < 1) return;
  if (diag->nxings == 1) {
    diag->nxings = 0;
    diag->is_link = true;
    return;
  }

  // diagram has at least two xings

  
}

void deleteXing (DiagramPtr diag, XingPtr xp) {
  if (superv) 
    fprintf (stderr, "deleting Xing %d from %d\n", xp->ID, which_reid_move);
  static int order = 1;
  if (fpkps) {
    if (which_reid_move == RI) {
      fprintf (fpkps, "txt fcolor #d00 .35\n");
      fprintf (fpkps, "txt %f %f 5 \"RI %d\"\n", xp->x, xp->y, order);
    }
    else if (which_reid_move == RII){
      fprintf (fpkps, "txt fcolor #0d0 .35\n");
      fprintf (fpkps, "txt %f %f 5 \"RII %d\"\n", xp->x, xp->y, order);
    }
    else if (which_reid_move == RRC) {
    }
    else
      panic_exit ("deleteXing(): impossible case");

    ++order;
  }

  // move the Xing about to be deleted to the end of the Xing array

  swapXings (diag, xp->location, diag->nxings - 1);

  if (xp == diag->first_xing)
    diag->first_xing = diag->xing [0];
  diag->nxings--;
  diag->removals++;
}

bool deleteNugatoryXing (DiagramPtr diag, XingPtr xp) {
  checkDiagram (diag, "deleteNugatoryXing A");
  int tstrand = 0;
  if (xp->xtop [0] == xp)
    tstrand = 1;
  int bstrand = 1 - tstrand;

  tie_ends (xp->xbot [bstrand], xp->sbot [bstrand], 
            xp->xtop [tstrand], xp->stop [tstrand]);

  deleteXing (diag, xp);

  checkDiagram (diag, "deleteNugatoryXing B");
  return true;
}

void recomboDiagram (DiagramPtr diag) {
  if (diag->nxings < 2) return;
  // choose a xing at random
  XingPtr xp = diag->xing [rand_integer (0, diag->nxings)];
  //XingPtr xp = diag->xing [0];
  if (xp->xtop [0] == xp || xp->xtop [1] == xp) return;

  tie_ends (xp->xbot [0], xp->sbot [0],
	    xp->xtop [1], xp->stop [1]);
  tie_ends (xp->xbot [1], xp->sbot [1],
	    xp->xtop [0], xp->stop [0]);
  which_reid_move = RRC;
  deleteXing (diag, xp);
}

void compute_emcode_abcd (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    if (xp->sign == 1) {
      xp->embot [0] = 'c';
      xp->embot [1] = 'b';
      xp->emtop [0] = 'a';
      xp->emtop [1] = 'd';
    }
    else {
      xp->embot [0] = 'd';
      xp->embot [1] = 'c';
      xp->emtop [0] = 'b';
      xp->emtop [1] = 'a';
    }
  }
}

void EwingMillettCode (FILE *fp, DiagramPtr diag) {
  if (!labelDiagram (diag)) panic_exit ("EwingMillettCode(): diagram couldn't be labelled");
  
  compute_emcode_abcd (diag);

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    fprintf (fp, "%d", i + 1);
    if (xp->sign > 0) {
      fprintf (fp, "+");
      fprintf (fp, "%d%c", xp->xtop [0]->location + 1, xp->xtop [0]->embot [xp->stop [0]]); 
      fprintf (fp, "%d%c", xp->xbot [1]->location + 1, xp->xbot [1]->emtop [xp->sbot [1]]); 
      fprintf (fp, "%d%c", xp->xbot [0]->location + 1, xp->xbot [0]->emtop [xp->sbot [0]]); 
      fprintf (fp, "%d%c", xp->xtop [1]->location + 1, xp->xtop [1]->embot [xp->stop [1]]); 
    }
    else {
      fprintf (fp, "-");
      fprintf (fp, "%d%c", xp->xtop [1]->location + 1, xp->xtop [1]->embot [xp->stop [1]]); 
      fprintf (fp, "%d%c", xp->xtop [0]->location + 1, xp->xtop [0]->embot [xp->stop [0]]); 
      fprintf (fp, "%d%c", xp->xbot [1]->location + 1, xp->xbot [1]->emtop [xp->sbot [1]]); 
      fprintf (fp, "%d%c", xp->xbot [0]->location + 1, xp->xbot [0]->emtop [xp->sbot [0]]); 
    }

    fprintf (fp, "\n");
    
    
  }
}

void checkDiagram (DiagramPtr diag, char *calledfrom) {
#ifdef TRYING_TO_FIND_BUG
  //fprintf (stderr, "checking diagram...\n"); fflush (stderr);

  
  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];

    if (xp->location != i) 
      panic_exit ("checkDiagram: Xing %d bad location", i);

    if (xp->xtop [0]->xbot [xp->stop [0]] != xp) panic_exit ("checkDiagram from %s: top 0, nxings = %d", calledfrom, diag->nxings);
    if (xp->xtop [1]->xbot [xp->stop [1]] != xp) panic_exit ("checkDiagram from %s: top 1, nxings = %d", calledfrom, diag->nxings);
    if (xp->xbot [0]->xtop [xp->sbot [0]] != xp) panic_exit ("checkDiagram from %s: bot 0, nxings = %d", calledfrom, diag->nxings);
    if (xp->xbot [1]->xtop [xp->sbot [1]] != xp) panic_exit ("checkDiagram from %s: bot 1, nxings = %d", calledfrom, diag->nxings);
    
  }
  labelDiagram (diag);
  labelDiagram2 (diag);
#endif
}

bool removeXingsR1 (DiagramPtr diag) {
  if (diag->nxings <= diag->min_nxings) return false;

  checkDiagram (diag, "removeXingsR1 A");

  // returns true if it's possible that more Xings might be removed
  if (diag->zok) {
    if (diag->nxings < 3) {
      if (fpkps) {
        while (diag->nxings) 
          deleteNugatoryXing (diag, diag->first_xing);
      }
      else
        diag->nxings = 0;
      return false;
    }
  }
  else if (diag->nxings < 2)
    return false;
  
  checkDiagram (diag, "removeXingsR1 B");


  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    if (xp->xtop [0] == xp || xp->xtop [1] == xp) {
      return deleteNugatoryXing (diag, xp);
    }
  }

  checkDiagram (diag, "removeXingsR1 C");

  return false;
}


bool removeXingsR2 (DiagramPtr diag) {
  if (diag->nxings <= diag->min_nxings) return false;

  checkDiagram (diag, "removeXingsR2 A");

  // returns true if it's possible that more Xings might be removed

  if (diag->nxings < 3)   // let removeXingsR1() handle it
    return false;

  for (int uu = 0; uu < diag->nxings; uu++) {
    if (isNugatory (diag->xing [uu]))
      panic_exit ("removeXingsR2(): diagram has nugatory crossings!");
  }

  checkDiagram (diag, "removeXingsR2 A");

  // this subroutine assumes it is not possible to remove any Xings via R1 moves!!

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    XingPtr xpo;

    // case strands are parallel
    if (xp->xtop [0] == xp->xtop [1] && xp->sign != xp->xtop [0]->sign) {
      xpo = xp->xtop [0];
      if (xpo->xbot [0] != xp) panic_exit ("removeXingsR2(): impossible case A1");
      if (xpo->xbot [1] != xp) panic_exit ("removeXingsR2(): impossible case A2");

      tie_ends (xp->xbot [0], xp->sbot [0], xp->xtop [0]->xtop [1], xp->xtop [0]->stop [1]);
      tie_ends (xp->xbot [1], xp->sbot [1], xp->xtop [0]->xtop [0], xp->xtop [0]->stop [0]);
      deleteXing (diag, xp->xtop [0]);
      deleteXing (diag, xp);

      checkDiagram (diag, "removeXingsR2 B");
      return true;
    }

    // case strands circulate about centre in a CW fashion
    if (xp->xtop [0] == xp->xbot [1] && xp->sign != xp->xtop [0]->sign) {
      xpo = xp->xtop [0];
      if (xpo->xtop [0] == xp && xpo->xbot [1] == xp) {

        tie_ends (xp->xbot [0], xp->sbot [0], xpo->xtop [1], xpo->stop [1]);
        tie_ends (xpo->xbot [0], xpo->sbot [0], xp->xtop [1], xp->stop [1]);
        deleteXing (diag, xpo);
        deleteXing (diag, xp);

        checkDiagram (diag, "removeXingsR2 C");
        return true;
      }
    }

    // case strands circulate about centre in a CCW fashion
    if (xp->xtop [1] == xp->xbot [0] && xp->sign != xp->xtop [1]->sign) {
      xpo = xp->xtop [1];
      if (xpo->xtop [1] == xp && xpo->xbot [0] == xp) {

        tie_ends (xp->xbot [1], xp->sbot [1], xp->xtop [1]->xtop [0], xp->xtop [1]->stop [0]);
        tie_ends (xp->xtop [1]->xbot [1], xp->xtop [1]->sbot [1], xp->xtop [0], xp->stop [0]);
        deleteXing (diag, xp->xtop [1]);
        deleteXing (diag, xp);

        checkDiagram (diag, "removeXingsR2 D");
        return true;
      }
    }
  }

  return false;
}

XingPtr isIsolatedClasp (XingPtr xp) {
  XingPtr xn;

  // NOTE: only a xing that is on the bottom end of a clasp will return a non-NULL value

  // case strands are parallel
  if (xp->xtop [0] == xp->xtop [1] && xp->sign == xp->xtop [0]->sign) {
    xn = xp->xtop [0];
    if (xn == xp) return (XingPtr) NULL;
    if (xp->xbot [0] == xp->xbot [1]) return (XingPtr) NULL;
    if (xn->xtop [0] == xn->xtop [1]) return (XingPtr) NULL;
    return xn;
  }

  // case strands circulate CW
  if (xp->xtop [0] == xp->xbot [1] && xp->sign == xp->xtop [0]->sign) {
    xn = xp->xtop [0];
    if (xn->xtop [0] == xp && xn->xbot [1] == xp) {
      if (xn == xp) return (XingPtr) NULL;
      if (xp->xbot [0] == xp->xtop [1]) return (XingPtr) NULL;
      if (xn->xbot [0] == xn->xtop [1]) return (XingPtr) NULL;
      return xn;
    }
  }

  // case strands circulate CCW
  if (xp->xtop [1] == xp->xbot [0] && xp->sign == xp->xtop [1]->sign) {
    xn = xp->xtop [1];
    if (xn->xtop [1] == xp && xn->xbot [0] == xp) {
      if (xn == xp) return (XingPtr) NULL;
      if (xp->xbot [1] == xp->xtop [0]) return (XingPtr) NULL;
      if (xn->xbot [1] == xn->xtop [0]) return (XingPtr) NULL;
      return xn;
    }
  }
  return (XingPtr) NULL;
}


bool isNugatory (XingPtr xp) {
  if (xp->xtop [0] == xp || xp->xtop [1] == xp) 
    return true;
  return false;
}

XingPtr isSuperXing (XingPtr xp) {
  if (xp->kind & XING_KIND_ISOCLASP || isNugatory (xp))
    return (XingPtr) NULL;

  if (xp->xtop [0] == xp->xtop [1] && xp->sign == xp->xtop [0]->sign)
    return xp;

  if (xp->xbot [0] == xp->xbot [1] && xp->sign == xp->xbot [0]->sign)
    return xp;

  if (xp->xbot [0] == xp->xtop [1] && xp->sign == xp->xbot [0]->sign)
    return xp;

  if (xp->xtop [0] == xp->xbot [1] && xp->sign == xp->xtop [0]->sign)
    return xp;

  return (XingPtr) NULL;
}

bool not_working_isRIIremovable (XingPtr xp) {   // may not be correct

  XingPtr xn;
  if (xp->xtop [0] == xp->xtop [1] && xp->sign != xp->xtop [0]->sign) {
    xn = xp->xtop [0];
    xp->xother = xn;
    if (xn == xp) return false;
    return true;
  }
  else if (xp->xtop [0] == xp->xbot [1] && xp->sign != xp->xtop [0]->sign) {
    xn = xp->xtop [0];
    xp->xother = xn;
    if (xn == xp) false;
    return true;
  }
  else if (xp->xtop [1] == xp->xbot [0] && xp->sign != xp->xtop [1]->sign) {
    xn = xp->xtop [1];
    xp->xother = xn;
    if (xn == xp) return false;
    return true;
  }
  return false;
}


#define TRUE  1
#define FALSE 0

int findIsolatedClasps (DiagramPtr diag) {
  XingPtr xp, xn;
  int kount = 0;
  for (int i = 0; i < diag->nxings; i++) {
    diag->xing [i]->mark = FALSE;
    diag->xing [i]->kind = XING_KIND_NORMAL;
    if (diag->xing [i]->location != i) 
      panic_exit ("findIsolatedClasps(): impossible case A!");
  }

  
  for (int i = 0; i < diag->nxings; i++) {
    xp = diag->xing [i]; 
    if (xn = isIsolatedClasp (xp)) {
      if (!xp->mark && !xn->mark) {
        //printf ("%d %d isolated clasp\n", xp->ID, xn->ID);  // leave in for now
        xp->mark = xn->mark = TRUE;
	xp->kind = xn->kind = XING_KIND_ISOCLASP;
        diag->scratch [2 * kount] = xp->location;
        diag->scratch [2 * kount + 1] = xn->location;
        kount++;
      }
      else if ((xp->mark && !xn->mark) || (!xp->mark && xn->mark))
        panic_exit ("findIsolatedClasps(): impossible case B!");

    }
  }
  return kount;
}

int findSuperXings (DiagramPtr diag, char which) {
  XingPtr xp;
  int kount = 0;
  int sign = 1;
  findIsolatedClasps (diag);

  if (which == FLIP_NEG_SUPERXING)
      sign = -1;

  for (int i = 0; i < diag->nxings; i++) {
    xp = diag->xing [i]; 
    if (isSuperXing (xp) && xp->sign == sign) {   // shouldn't this be isSuperXing()?
      diag->scratch [kount] = xp->location;
      kount++;
    }
  }
  return kount;
}


int findPosNegXings (DiagramPtr diag, char which) {
  XingPtr xp;
  int kount = 0;
  int sign = 1;

  if (which == FLIP_NEG_XING)
      sign = -1;

  for (int i = 0; i < diag->nxings; i++) {
    xp = diag->xing [i]; 
    if (xp->sign == sign) {
      diag->scratch [kount] = xp->location;
      kount++;
    }
  }
  return kount;
}


int findNonIsolatedClasps (DiagramPtr diag) {
  XingPtr xp;
  int kount = 0;
  
  for (int i = 0; i < diag->nxings; i++) {
    xp = diag->xing [i]; 
    if (!isIsolatedClasp (xp)) {
      diag->scratch [kount] = xp->location;
      kount++;
    }
  }
  return kount;
}



int printIsolatedClasps (FILE *fp, DiagramPtr diag) {
  //fprintf (fp, "\n\n");
  int kount = findIsolatedClasps (diag);
  fprintf (fp, "%d isolated clasp", kount);
  if (kount != 1) printf ("s");
  if (kount > 0) {
    for (int i = 0; i < kount; i++) {
      fprintf (fp, " %d/%d", diag->xing [diag->scratch [2 * i]]->ID, diag->xing [diag->scratch [2 * i + 1]]->ID);
    }
  }
  fprintf (fp, "\n");
  return kount;
}

void output_scenery (FILE *fpscenery, FILE *fpscript, XingPtr xp1, XingPtr xp2, int label) {
  fprintf (fpscript, "txt %f %f 5 %d\n", xp1->x, xp1->y, label);
  fprintf (fpscript, "txt %f %f 5 %d\n", xp2->x, xp2->y, label);
}

void classifyXings (DiagramPtr diag, int &niso, int &nnug, int &nsup) {
  niso = 2 * findIsolatedClasps (diag);
  nnug = nsup = 0;
  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];   // for convenience
    //xp->kind = XING_KIND_NORMAL;

    // note: it is an error if a xing is classified as more than one type
    // xp->kind is ORed together to allow detection of this error

    //if (isIsolatedClasp (xp))
    //xp->kind |= XING_KIND_ISOCLASP;
    if (isNugatory (xp)) {
      nnug++;
      xp->kind |= XING_KIND_NUGATORY;
    }
    if (isSuperXing (xp)) {
      nsup++;
      xp->kind |= XING_KIND_SUPERXING;
    }
  }
}

char xing_kind_name [3][12] = {
  "nugatory",
  "isoclasp",
  "superxing"};

void output_KnotPlot_script_and_scenery (char *basename, DiagramPtr diag) {
  char *script = (char *) calloc (strlen (basename) + 8, sizeof (char));
  sprintf (script, "%s.kps", basename);
  FILE *fpscript = open_file (script, "w");

  int iso_kount = findIsolatedClasps (diag);
  //KnotGroupRepresentation (diag);
  
  int iso_kount2 = 0;
  int super_kount = 0;
  int nugatory_kount = 0;
  int normal_kount = 0;

  int nsup, nnug, niso;
  classifyXings (diag, niso, nnug, nsup);

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    if (xp->kind == XING_KIND_NORMAL)
      ++normal_kount;
    if (xp->kind & XING_KIND_ISOCLASP)
      ++iso_kount2;
    if (xp->kind & XING_KIND_NUGATORY)
      ++nugatory_kount;
    if (xp->kind & XING_KIND_SUPERXING)
      ++super_kount;
  }

  bool ok = true;

  if (nsup != super_kount || niso != iso_kount2 || nnug != nugatory_kount) {
    fprintf (stderr, "Warning: classification error %d %d  %d %d  %d %d\n",
	     nsup, super_kount, niso, iso_kount2, nnug, nugatory_kount);
    ok = false;
  }

  if (normal_kount + iso_kount2 + nugatory_kount + super_kount != diag->nxings) {
    fprintf (stderr, "Warning: xing count mismatch: %d + %d + %d + %d != %d\n",
	     normal_kount, iso_kount2, nugatory_kount, super_kount, diag->nxings);
    ok = false;
  }
  
  if (2 * iso_kount != iso_kount2) {
    fprintf (stderr, "Warning: clasp count mismatch: 2 * %d != %d\n", 
	     iso_kount, iso_kount2);
    ok = false;
  }
  
 			   

  fprintf (fpscript, "ortho\n%%txt init\ntxt delete group 9\ntxt delete group 7\n;txt delete group 11;txt font 4\ntxt centre\n");
  fprintf (fpscript, "txt push;txt group 11;txt c -20 \"crossing classification\n");
  fprintf (fpscript, "txt font 3\ntxt c -40 \"xinger version %s %s\"\ntxt pop\ntxt on\n\n", __DATE__, __TIME__);
  fprintf (fpscript, "scenery init\nscenery dscale 1.5\n3dop ostep 3\n3dop oscale .5\nshow orient\n\n");

  if (!ok) {
    fprintf (fpscript, "txt font 4\ntxt c c failed\n");
    fclose (fpscript);
    panic_exit ("failed for some reason");
  }


  int kind = 1;
  double z = 5.0;
  double zinc = 1.08;

  fprintf (fpscript, "silent=t\n");
  for (int i = 0; i < 3; i++)  {
    for (int x = 0; x < diag->nxings; x++) {
      char extra;
      if (diag->xing [x]->kind & kind) {
	if (kind == XING_KIND_SUPERXING && diag->xing [x]->sign == 1) 
	  extra = '+';
	else if (kind == XING_KIND_SUPERXING && diag->xing [x]->sign == -1) 
	  extra = '-';
	else
	  extra = ' ';

	fprintf (fpscript, "scenery object %s%c %f %f %f\n", xing_kind_name [i], extra, diag->xing [x]->x, diag->xing [x]->y, z);
      }
    }
    kind *= 2;
    z += zinc;
  }
  fprintf (fpscript, "silent=f\n");
  real xoff, yoff;
  fprintf (fpscript, "\ntxt push\ntxt group 9\ntxt cent\ntxt font 4\n");
  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    fprintf (fpscript, "\n%% xing #%d, sign = %d, angle = %f, xtop [1] = %d, xtop [0] = %d, xbot [1] = %d, xbot [0] = %d\n", 
	     i + 1, xp->sign, xp->angle * RAD_TO_DEG,
	     xp->xtop [1]->location + 1, xp->xtop [0]->location + 1, xp->xbot [1]->location + 1, xp->xbot [0]->location + 1);
    fprintf (fpscript, "txt %f %f %f %d\n",  diag->xing [i]->x, diag->xing [i]->y, z, i + 1);
    // xp->angle is 45 degree CCW rotation of overcrossing angle
    
    xoff = cos (xp->angle);
    yoff = sin (xp->angle);

    if (xp->sign == 1)
      fprintf (fpscript, "txt %f %f %f +\n",  diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z);
    else
      fprintf (fpscript, "txt %f %f %f -\n",  diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z);
  }
  fprintf (fpscript, "\ntxt pop\n");

  fprintf (fpscript, "\ntxt push\ntxt group 7\ntxt cent\ntxt font 3\n");

  compute_emcode_abcd (diag);

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    // xp->angle is 45 degree CCW rotation of overcrossing angle
    real angle = xp->angle;
    if (xp->sign == 1) 
      angle += 45.0 * DEG_TO_RAD;
    else 
      angle -= 45.0 * DEG_TO_RAD;

    xoff = cos (angle);
    yoff = sin (angle);
    fprintf (fpscript, "txt %f %f %f %c\n",  
	     diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z,
	     xp->emtop [1]);
    angle -= 90.0 * DEG_TO_RAD;
    xoff = cos (angle); yoff = sin (angle);
   
    fprintf (fpscript, "txt %f %f %f %c\n",  
	     diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z,
	     xp->emtop [0]);
    angle -= 90.0 * DEG_TO_RAD;
    xoff = cos (angle); yoff = sin (angle);
   
    fprintf (fpscript, "txt %f %f %f %c\n",  
	     diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z,
	     xp->embot [1]);
    angle -= 90.0 * DEG_TO_RAD;
    xoff = cos (angle); yoff = sin (angle);
   
    fprintf (fpscript, "txt %f %f %f %c\n",  
	     diag->xing [i]->x + xoff, diag->xing [i]->y + yoff, z,
	     xp->embot [0]);
  }
  fprintf (fpscript, "\ntxt pop\n");


  fprintf (fpscript, "\n%%scenery on\n");
  fclose (fpscript);
}

bool flipXing (DiagramPtr diag) {
  if (diag->nxings < 1) return false;
  XingPtr xp = diag->xing [rand_integer (0, diag->nxings)];
  xp->sign *= -1;
  diag->xing_flips++;
  if (fp_ray)
    fprintf (fp_ray, "n %f %f %f\n", xp->x, xp->y, xp->z);
  return true;
}

bool flip_first = false;

bool flipIsolatedClasp (DiagramPtr diag) {
  int kount = findIsolatedClasps (diag);
  diag->iso_clasps += kount;

  if (kount < 1) {
    if (flip_first)
      return false;
    else
      return flipXing (diag);
  }

  int flipA = 2 * rand_integer (0, kount);

  if (flip_first) 
    flipA = 0;

  
  int a = flipA;
  int b = flipA + 1;
  if (rand_uniform () < 0.5) 
    flipA++;

  XingPtr xp = diag->xing [diag->scratch [flipA]];
  xp->sign *= -1;
  if (fp_ray)
    fprintf (fp_ray, "i %f %f %f\n", xp->x, xp->y, xp->z);

  //printf ("Xing #%d got flipped\n", xp->ID);
  diag->iso_clasp_flips++;
  return true;
}


bool flipSuperXing (DiagramPtr diag, char which) {
  int kount = findSuperXings (diag, which);
  if (kount < 1) return false;

  int flipA = rand_integer (0, kount);
  
  XingPtr xp = diag->xing [diag->scratch [flipA]];
  xp->sign *= -1;

  return true;
}


bool flipPosNegXing (DiagramPtr diag, char which) {
  int kount = findPosNegXings (diag, which);
  if (kount < 1) return false;

  int flipA = rand_integer (0, kount);
  
  XingPtr xp = diag->xing [diag->scratch [flipA]];
  xp->sign *= -1;

  return true;
}



bool superFlip (DiagramPtr diag) {
  int kount = findNonIsolatedClasps (diag);
  if (kount < 1) return false;
  int flipA = rand_integer (0, kount);
  
  XingPtr xp = diag->xing [diag->scratch [flipA]];
  xp->sign *= -1;

  return true;
}

/*
what should the program do if there are no isolated clasps and p=0?
*/

void reflectZ (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) 
    diag->xing [i]->sign *= -1;
}


void removeXings (DiagramPtr diag, bool R2) {
  int nxings_before;
  do {
    nxings_before = diag->nxings;

    checkDiagram (diag, "removeXings A");

    // always remove via R1 (for now)
    which_reid_move = RI;
    while (removeXingsR1 (diag));  // empty loop

    checkDiagram (diag, "removeXings B");

    which_reid_move = RII;
    if (R2)
      removeXingsR2 (diag);
    checkDiagram (diag, "removeXings C");
  } while (nxings_before != diag->nxings);

    checkDiagram (diag, "removeXings D");

}

void writeNumXings (FILE *fp, DiagramPtr diag) {
  fprintf (fp, "%d\n", diag->nxings);
} 

void writeNumClasps (FILE *fp, DiagramPtr diag) {
  int kount = findIsolatedClasps (diag);
  fprintf (fp, "%d\n", kount);
} 

void flipAllXings (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) {
    if (rand_uniform () < 0.5) 
      diag->xing [i]->sign = 1;
    else
      diag->xing [i]->sign = -1;
  }
}

int gumby = 0;

void output_colours (DiagramPtr diag) {
  FILE *fp = fopen ("c.kps", "w");
  fprintf (fp, "scenery init\n");
  for (int i = 0; i < diag->nxings; i++) {
    if (diag->xing [i]->colour) {
      fprintf (fp, "scenery object %s %f %f %f\n", diag->xing [i]->colour,
	       diag->xing [i]->x, diag->xing [i]->y, 0.0);
      if (!strcmp ("white", diag->xing [i]->colour)) {
	fprintf (fp, "txt init\ntxt on\ntxt font 4\ntxt cent\n");
	XingPtr xp = diag->xing [i];
	fprintf (fp, "txt %f %f 4.0 top0\n", xp->xtop [0]->x, xp->xtop [0]->y);
	fprintf (fp, "txt %f %f 4.0 top1\n", xp->xtop [1]->x, xp->xtop [1]->y);
	fprintf (fp, "txt %f %f 4.0 bot0\n", xp->xbot [0]->x, xp->xbot [0]->y);
	fprintf (fp, "txt %f %f 4.0 bot1\n", xp->xbot [1]->x, xp->xbot [1]->y);
      }
    }
  }
  fprintf (fp, "scenery on\n");
  fprintf (fp, "txt c -20 g%d\n", gumby);
  fclose (fp);

  return;
  XingPtr xp, rxp1, lxp1, rxp0, lxp0;	
  int dum;
  for (int i = 0; i < diag->nxings; i++) {
    xp = diag->xing [i];
    right (xp, 0, rxp0, dum);
    right (xp, 1, rxp1, dum);
    left (xp, 0, lxp0, dum);
    left (xp, 1, lxp1, dum);

    printf ("%d: %d %d %d %d  %d %d  %d %d\n",
	    xp->gauss,
	    xp->xbot [0]->gauss,
	    xp->xbot [1]->gauss,
	    xp->xtop [0]->gauss,
	    xp->xtop [1]->gauss,
	    lxp0->gauss, rxp0->gauss,
	    lxp1->gauss, rxp1->gauss
	    );
  }
}


bool colour_it = false;

void colourXing (XingPtr xp, char *colour) {
  xp->colour = (char *) calloc (strlen (colour) + 2, sizeof (char));
  strcpy (xp->colour, colour);
}

bool isSuperNugatory (XingPtr xp, int strand) {
  XingPtr xp_start = xp;
  int strand_start = strand;
  fore (xp_start, strand_start, xp, strand);
  if (xp == xp_start)  // xing is ordinary nugatory
    return false;

  int ab = abovebelow (xp, strand);

  do {
    if (abovebelow (xp, strand) != ab) return false;
    fore (xp, strand, xp, strand);
  } while (xp != xp_start);

  if (strand == strand_start) 
    panic_exit ("isSuperNugatory(): impossible case!");

  return true;
}

int entangling (DiagramPtr diag) {
  if (diag->nxings < 2)
    return 0;

  int kount = 0;
  XingPtr xp = diag->first_xing;
  int strand = xp->start;
  
  do {
    XingPtr fxp;
    int fstrand;
    fore (xp, strand, fxp, fstrand);
    if (abovebelow (xp, strand) != abovebelow (fxp, fstrand)) kount++;
    xp = fxp;
    strand = fstrand;
  } while (!(xp == diag->first_xing && strand == xp->start));

  return kount;
}

int longestBridge (DiagramPtr diag) {
  if (diag->nxings < 2)
    return 1;

  int max = 0;

  return max;
}

bool detectR3R2simps (XingPtr xp, int strand) {
  //colour_it = true;
  XingPtr fxp, rxp, lxp, rfxp, lfxp, rbxp, lbxp, bxp;
  int fstrand, rstrand, lstrand, rfstrand, lfstrand, rbstrand, lbstrand, bstrand;

  fore (xp, strand, fxp, fstrand);
  back (xp, strand, bxp, bstrand);
  //

  right (xp, strand, rxp, rstrand);
  left (xp, strand, lxp, lstrand);

  right (fxp, fstrand, rfxp, rfstrand);
  left (fxp, fstrand, lfxp, lfstrand);

  // first check to see if any R3 R1 removals are possible

  if ((lxp == lfxp && rxp == fxp && abovebelow (xp, strand) != abovebelow (lxp, lstrand)) ||
      (rxp == rfxp && lxp == fxp && abovebelow (xp, strand) != abovebelow (rxp, rstrand))) 
    return true;

  if (isSuperNugatory (xp, strand)) 
    return true;

  static bool x = false;
  if (x) {
    x = false;
    
    printf ("strand %d, fstrand %d, lstrand %d, rstrand %d\n",
	    strand, fstrand, lstrand, rstrand);
    colourXing (xp, "white");
    colourXing (fxp, "blue");
    colourXing (rxp, "red");
    colourXing (lxp, "green");
    colourXing (bxp, "yellow");
  }

  

  if (rxp != rfxp || lxp != lfxp) return false;

  if (abovebelow (rxp, rstrand) != abovebelow (lxp, lstrand)) return false;
  if (abovebelow (xp, strand) != abovebelow (fxp, fstrand) &&
      abovebelow (xp, strand) == abovebelow (lxp, lstrand)) return false;
  

  /*
  right (bxp, bstrand, rbxp, rbstrand);
  left (bxp, bstrand, lbxp, lbstrand);

  if (rxp != rbxp || lxp != lbxp) return false;
  */
  if (rstrand == rfstrand || lstrand == lfstrand) 
    panic_exit ("detectrR3R2simps(): impossible case!");

  if (colour_it) {
    colourXing (xp, "white");
    colourXing (fxp, "yellow");
    colourXing (rxp, "red");
    colourXing (lxp, "green");
  }

  return true;
}



void detectR3R2simps (DiagramPtr diag) {
  // note this doesn't simplify the knot and maintain knot type,
  // it just sets number of xings to 0 so that the diagrams is 
  // not considered in the search for new monsters

  // only works for proper knots
  XingPtr xp = diag->first_xing;
  //xp->start = gumby;
  int strand = xp->start;

  // printf ("starting at %d\n", strand);
  do {
    if (detectR3R2simps (xp, strand)) {
      if (!colour_it) 
	diag->nxings = 0;
      return;
    }
    XingPtr xpn;
    xpn = xp->xtop [strand];
    strand = xp->stop [strand];
    xp = xpn;
  } while (!(xp == diag->first_xing && strand == xp->start));

}

void printDebugInfo (DiagramPtr diag, FILE *fp) {
  if (!labelDiagram (diag)) panic_exit ("printDebugInfo(): diagram couldn't be labelled");

  for (int i = 0; i < diag->nxings; i++) {
    XingPtr xp = diag->xing [i];
    fprintf (fp, "%d: ID %d g %d t0 %d:%d t1 %d:%d b0 %d:%d b1 %d:%d\n",
	     i,
	     xp->ID, xp->gauss,
	     xp->xtop [0]->ID, xp->stop [0],
	     xp->xtop [1]->ID, xp->stop [1],
	     xp->xbot [0]->ID, xp->sbot [0],
	     xp->xbot [1]->ID, xp->sbot [1]);
  }
}
	     
    

void executeProgram (DiagramPtr diag, char *program, FILE *fp) {
  char *command = program;
  
  int save_min = diag->min_nxings;

  while (*command) {
    switch (*command) {
    case PRINT_DEBUG_INFO:
      printDebugInfo (diag, fp);
      break;
    case PRINT_ENTANGLING:
      fprintf (fp, "%d ", entangling (diag));
      break;
    case PRINT_BRIDGE_INFO:
      fprintf (fp, "%d ", longestBridge (diag));
      break;
    case DETECT_R3_R2_SIMPS:  // option should always be last operator in program
      detectR3R2simps (diag); // also, should be used with the -znb option
                              // NOTE: this option DOES NOT maintain knot type, it simply
                              // detects if there are R3R2simps and if there are it
                              // sets the number of crossings to zero.
      break;
    case STANDARD_LABELLING:
      findStandardLabelling (diag);
      break;
    case WRITE_NUM_XINGS:
      writeNumXings (stdout, diag);
      break;
    case WRITE_NUM_CLASPS:
      writeNumClasps (stdout, diag);
      break;
    case RECOMBO:
      support_links = true;
      recomboDiagram (diag);
      break;
    case FLIP_ALL_XINGS:
      flipAllXings (diag);
      break;
    case CHECK_DIAGRAM:
      checkDiagram (diag, "executeProgram A");
      break;
    case REIDREM_via_R1R2_0:
      diag->min_nxings = 0;   // fall thru
    case REIDREM_via_R1R2:
      removeXings (diag, true);
      diag->min_nxings = save_min;
      break;
    case REIDREM_via_R1:
      checkDiagram (diag, "executeProgram D");
      removeXings (diag, false);
      checkDiagram (diag, "executeProgram E");
      break;
    case FLIP_XING:
      flipXing (diag);
      break;
    case SUPER_FLIP:
      superFlip (diag);
      break;
    case FLIP_POS_SUPERXING:
    case FLIP_NEG_SUPERXING:
      if (rand_uniform () < diag->p)
	flipSuperXing (diag, *command);
      else 
	flipXing (diag);
      break;
    case FLIP_POS_XING:
    case FLIP_NEG_XING:
      flipPosNegXing (diag, *command);
      break;
    case FLIP_FIRST_ISOCLASP:
      flip_first = true;  
      diag->p = 1.08;    // fall thru (i.e., no break)
    case CLASP_SP:
      checkDiagram (diag, "executeProgram F");
      if (rand_uniform () < diag->p)
        flipIsolatedClasp (diag);
      else 
        flipXing (diag);
      checkDiagram (diag, "executeProgram G");
      break;
    case PRINT_ISO_CLASPS:
      printIsolatedClasps (stdout, diag);
      break;
    case REVERSE_ORIENTATION:
      panic_exit ("reverse orientation not yet implemented");
      break;
    case REFLECT_IN_Z:
      reflectZ (diag);
      break;
    default:
      panic_exit ("unknown operator `%c'\n", *command);
    }
    command++;
  }
}

void standardLabelling (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) {
    for (int s = 0; s < 2; s++) {
      diag->first_xing = diag->xing [i];
      diag->first_xing->start = s;
      GaussCode (stderr, diag);
    }
  }
}

void findStandardLabelling (DiagramPtr diag) {
  if (diag->nxings < 2) return;
  strcpy (diag->cscratch, "zzzz");  
  int ssave = 108; int isave = 102;  // these values always get overwritten

  for (int i = 0; i < diag->nxings; i++) {
    for (int s = 0; s < 2; s++) {
      diag->first_xing = diag->xing [i];
      diag->first_xing->start = s;
      GaussCode ((FILE *) NULL, diag);
      if (strcmp (diag->cscratch, diag->gauss_code) > 0) {
	strcpy (diag->cscratch, diag->gauss_code);
	isave = i;
	ssave = s;
      }
    }
  }
  diag->first_xing = diag->xing [isave];
  diag->first_xing->start = ssave;
}

bool wag2 = false;
void uuu (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) {
    if (diag->xing [i]->sign == 1) 
      printf ("1");
    else
      printf ("0");
  }
  printf ("\n");
}

void AllGaussCodes (FILE *fp, DiagramPtr diag) {
  if (diag->nxings > 30) 
    panic_exit ("knot has too many crossings!  Maximum is 30");

  int N = (int) (pow (2.0, diag->nxings) + 0.5);
  if (!wag2) N /= 2;

  for (int i = 0; i < N; i++) {
    int j = i;
    for (int x = 0; x < diag->nxings; x++) {
      if (0x01 & j) 
	diag->xing [x]->sign = 1;
      else 
	diag->xing [x]->sign = -1;
      j = j >> 1;
    }
    //uuu (diag);
    GaussCode (fp, diag);
  }
}

void NGaussCodes (FILE *fp, DiagramPtr diag, int N) {
  for (int i = 0; i < N; i++) {
    for (int x = 0; x < diag->nxings; x++) {
      if (rand_uniform () < 0.5) 
	diag->xing [x]->sign = 1;
      else 
	diag->xing [x]->sign = -1;
    }
    GaussCode (fp, diag);
  }
}

#define NO_REGION 0     // meaning that no region has been assigned

void markRegion_orig (int region, XingPtr xpstart, int c) {
  XingPtr xp = xpstart;

  // this routine will fail in some cases with nugatory crossings (see below)

  do {
    if (xp->region [c] != NO_REGION) 
      panic_exit ("no region panic");
    xp->region [c] = region;
    // find next xing and next c
    switch (c) {
    case 0:
      if (xp->sbot [0] == 0) 
	c = 1;
      else 
	c = 2;
      xp = xp->xbot [0];
      break;
    case 1:
      if (xp->sbot [1] == 0)
	c = 1;
      else 
	c = 2;
      xp = xp->xbot [1];
      break;
    case 2:
      if (xp->stop [0] == 0) 
	c = 3;
      else
	c = 0;
      xp = xp->xtop [0];
      break;
    case 3:
      if (xp->stop [1] == 0)
	c = 3;
      else
	c = 0;
      xp = xp->xtop [1];
      break;
    }
  } while (xp != xpstart);
}

void markRegion (int region, DiagramPtr diag, XingPtr xpstart, int c) {
  XingPtr xp = xpstart;

  // this function does an excessive amount of work in most cases, but
  // nugatory crossings mess things up (see above)

  int kount = 2 * diag->nxings;
  do {
    if (!(xp->region [c] == NO_REGION || xp->region [c] == region)) 
      continue;
    xp->region [c] = region;
    // find next xing and next c
    switch (c) {
    case 0:
      if (xp->sbot [0] == 0) 
	c = 1;
      else 
	c = 2;
      xp = xp->xbot [0];
      break;
    case 1:
      if (xp->sbot [1] == 0)
	c = 1;
      else 
	c = 2;
      xp = xp->xbot [1];
      break;
    case 2:
      if (xp->stop [0] == 0) 
	c = 3;
      else
	c = 0;
      xp = xp->xtop [0];
      break;
    case 3:
      if (xp->stop [1] == 0)
	c = 3;
      else
	c = 0;
      xp = xp->xtop [1];
      break;
    }
  } while (kount--);
}


void KnotGroupRepresentation (DiagramPtr diag) {
  for (int i = 0; i < diag->nxings; i++) {
    for (int c = 0; c < 4; c++) {
      diag->xing [i]->region [c] = NO_REGION;
    }
  }
  int region = 1;   // regions are numbered starting at 1
  
  for (int i = 0; i < diag->nxings; i++) {
    for (int c = 0; c < 4; c++) {
      if (diag->xing [i]->region [c] == NO_REGION) {
	markRegion (region++, diag, diag->xing [i], c);
      }
    }
  }
  
  diag->nregions = region - 1;
}

void KnotGroup (FILE *fp, DiagramPtr diag) {
  static int num = 0;
  if (diag->nxings < 3 && !diag->is_link) {
    fprintf (fp, "knot %d\n%8d%8d\n%6d%6d%6d%6d\nEND RECORD\n", num++, 3, 1, 3, -1, 2, -1);
    return;
  }
  KnotGroupRepresentation (diag);
  fprintf (fp, "knot %d\n%8d%8d\n", num++, diag->nregions, diag->nxings);
  for (int i = 0; i < diag->nxings; i++) {
    int sign = 1;
    int coff = 0;
    if (diag->xing [i]->sign == 1)
      coff = 3;
    for (int c = 0; c < 4; c++, sign *= -1) {
      fprintf (fp, "%6d", sign * diag->xing [i]->region [(c + coff) % 4]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "END RECORD\n");
}
