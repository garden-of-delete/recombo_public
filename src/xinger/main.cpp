#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#ifdef _WIN32   // disable annoying warnings
#pragma warning (disable: 4996)
#pragma warning (disable: 4101)
#pragma warning (disable: 4267)
#endif

#ifdef WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "xinger.h"
#include "random.h"
#include "util.h"

#define UNKNOWN    0
#define THREE_FILE 1
#define BINARY     2
#define GAUSS      3
#define IGNORE     4
#define KNOTPLOT_X 5
#define KNOTPLOT_S 6  // write KnotPlot script and KnotPlot scenery
#define JHOMFLY    7
#define ALLGAUSS   8
#define KGROUP     9
#define NUM_XINGS 10
#define NGAUSS    11
#define EMCODE    12

typedef struct {
  int dowker [2];
  int above_below [2];
  int sign [2];
} ExtendedGaussCode;


ExtendedGaussCode *egauss;

#define DEFAULT_MAX_XINGS 600

int input_format = UNKNOWN;
int output_format = GAUSS;

void printem (FILE *);

FILE *fpin = (FILE *) NULL;
FILE *fpout = (FILE *) NULL;
FILE *fpic = (FILE *) NULL, *fpid = (FILE *) NULL, *fpio = (FILE *) NULL;  // input file pointers for crossings, dowker and orientation, resp.
FILE *fpoc = (FILE *) NULL, *fpod = (FILE *) NULL, *fpoo = (FILE *) NULL;  // output file pointers for crossings, dowker and orientation, resp.

bool warn_deprecated = false;

DiagramPtr diag = (DiagramPtr) NULL;
int max_xings = DEFAULT_MAX_XINGS;
int *gauss;
int *sign;


Coord *coords;

int nxings;
int flags = 0;

#define FLAG_R1   1
#define FLAG_R2   2
#define FLAG_RREM (FLAG_R1 | FLAG_R2)

bool get_knot_three_file () {
  if (fscanf (fpic, "%d", &nxings) != 1) return false;
  if (nxings >= max_xings) panic_exit ("too many crossings (%d) found", nxings);
  int n2 = 2 * nxings;
  for (int i = 1; i <= n2; i++) 
    gauss [i] = 0;
  for (int i = 1; i <= n2; i += 2) {
    if (fscanf (fpid, "%d", gauss + i) != 1) panic_exit ("end of file before end of data on dowker");
    if (fscanf (fpio, "%d", sign  + i) != 1) panic_exit ("end of file before end of data on orient");
    int absgi = abs (gauss [i]);
    if (absgi % 2 != 0) panic_exit ("non-even dowker code found");
    if (absgi < 1 || absgi > n2) panic_exit ("out of bounds on dowker code");
    if (gauss [i] > 0)
      gauss [absgi] = -i;
    else
      gauss [absgi] = i;
    sign [absgi] = sign [i];
  }
  for (int i = 1; i <= n2; i++) 
    if (gauss [i] == 0) panic_exit ("not a valid dowker code");

  makeDiagramFromGauss (diag, gauss, sign, nxings);
  return true;
}

bool get_knot_KnotPlot_X () {
  static char *inp = (char *) 0;
  if (!inp) 
    inp = (char *) calloc (108, sizeof (char));

  // ignore first line (unless nothing is returned, in which case end of file has been reached)
  if (!fgets (inp, 102, fpin)) return false;
  //printf (">>%s<<\n", inp);
  if (!fgets (inp, 102, fpin)) panic_exit ("end of file before end of data");
  //printf (">>%s<<\n", inp);
  if (strncmp ("Found", inp, 5)) panic_exit ("not a KnotPlot-X file");
  sscanf (inp, "Found %d crossings", &nxings);
  if (nxings < 0 || nxings >= max_xings)
    panic_exit ("inappropriate number of crossings (%d) found!", nxings);

  int n2 = 2 * nxings;
  for (int i = 1; i <= n2; i++) 
    gauss [i] = 0;
  for (int i = 1; i <= n2; i += 2) {
    int ignore1, ignore2;
    double dignore1, dignore2, x, y, z, angle;

    if (!fgets (inp, 102, fpin)) panic_exit ("end of file before end of data");
    if (sscanf (inp, "%d%d%d%d%lf%lf%lf%lf%lf%lf)", 
		gauss + i, sign + i,
		&ignore1, &ignore2, 
		&dignore1, &dignore2, 
		&x, &y, &z, &angle) != 10) panic_exit ("unacceptable number of data items on input line");

    coords [i].x = x;
    coords [i].y = y;
    coords [i].z = z;
    coords [i].angle = angle;
    int absgi = abs (gauss [i]);
    if (absgi % 2 != 0) panic_exit ("non-even dowker code found");
    if (absgi < 1 || absgi > n2) panic_exit ("out of bounds on dowker code");
    if (gauss [i] > 0)
      gauss [absgi] = -i;
    else
      gauss [absgi] = i;
    sign [absgi] = sign [i];

    // skip a line
    fgets (inp, 102, fpin);
  }
  for (int i = 1; i <= n2; i++) 
    if (gauss [i] == 0) panic_exit ("not a valid dowker code");

  makeDiagramFromGauss (diag, gauss, sign, coords, nxings);
  return true;

}

void put_knot_three_file () {
  static bool first_call = true;
  if (first_call)
    first_call = false;
  else
    fprintf (fpoc, " ");

  DTcode (fpoc, fpod, fpoo, diag);
}

bool get_knot_binary () {
  unsigned short nxings_us;
  if (!read_unsigned_short (&nxings_us, fpin)) return false;
  nxings = (int) nxings_us;
  if (nxings >= max_xings)
    panic_exit ("too many crossings (%d) found", nxings);
  BitFile bf;
  bf.ch = (char) 0;
  bf.index = 0;
  bf.fp = fpin;
  int max_label = 4;
  int num_bits_per_xing = 1;
  while (2 * nxings > max_label){
    num_bits_per_xing++;
    max_label *= 2;
  }

  int n2 = 2 * nxings;
  for (int i = 1; i <= n2; i++) 
    gauss [i] = 0;

  // in the binary file, a num_bits_per_xing-bit binary number is used
  // to represent the DT code label, for example, the binary numbers
  // 000 001 010 011 100 ... represent the decimal DT code labels
  // 2   4   6   8   10  ... 

  // a repeated label indicates the start of a new link component

  for (int i = 1; i <= n2; i += 2) {
    sign [i] = read_bit (&bf) * 2 - 1;
    int dt_sign = read_bit (&bf) * 2 - 1;
    int absgi = 0;   
    int last_absgi;
    for (int b = 0; b < num_bits_per_xing; b++) {
      absgi |= read_bit (&bf) << b;
    }
    absgi = (absgi + 1) * 2;

    // if the DT label is the same as the previous one (not possible in a DT code),
    // start a new link component, have to do a i -= 2 at some point
    if (i > 1) {
      if (last_absgi == absgi)
        panic_exit ("reading binary files for links not yet implemented");
    }
    last_absgi = absgi;

    gauss [i] = absgi * dt_sign;

    if (gauss [i] > 0)
      gauss [absgi] = -i;
    else
      gauss [absgi] = i;
    sign [absgi] = sign [i];
  }

  for (int i = 1; i <= n2; i++) 
    if (gauss [i] == 0) panic_exit ("not a valid dowker code");

  int num_bits_read = nxings * (2 + num_bits_per_xing);
  // it's not really necessary to READ the pad bits (because of
  // the way read_bit() works), but it's a good check anyways.
  while (num_bits_read % 8 != 0) {
    if (read_bit (&bf) != 1)
      panic_exit ("pad bit not 1 as expected");
    num_bits_read++;
  }

  makeDiagramFromGauss (diag, gauss, sign, nxings);
  return true;
}

void put_knot_kgroup () {
  KnotGroup (fpout, diag);
}

void put_knot_binary () {
  DTcode (fpout, diag);
}

#define GAUSS_ab     2
#define GAUSS_DIGIT  3
#define GAUSS_SIGN   4

#define GAUSS_ABOVE  1
#define GAUSS_BELOW -1

void printem (FILE *fp) {
  int n2 = 2 * nxings;

  for (int i = 1; i <= n2; i++) 
    fprintf (fp, "%4d ", i);
  fprintf (fp, "\n");

  for (int i = 1; i <= n2; i++) 
    fprintf (fp, "%4d ", gauss [i]);
  fprintf (fp, "\n");

  for (int i = 1; i <= n2; i++) 
    fprintf (fp, "%4d ", sign [i]);
  fprintf (fp, "\n");
}

void add_gauss_crossing (int above_below, int label, int dowker_label, int sign2) {
  int index;
  if (egauss [label].dowker [0] == 0) 
    index = 0;
  else 
    index = 1;
  egauss [label].dowker [index] = dowker_label;
  egauss [label].above_below [index] = above_below;
  egauss [label].sign [index] = sign2;
}

void extendedGCtoDC (int max_label) {
  int index = 0;
  for (int j = 1; j <= nxings; j++, index++) {
    while (egauss [index].dowker [0] == 0) {
      if (index == max_xings || index > max_label) 
        panic_exit ("something horrible happened!");
      index++;
    }

    if (egauss [index].sign [0] != egauss [index].sign [1])
      panic_exit ("sign mismatch for label %d", index);
    if (egauss [index].above_below [0] == egauss [index].above_below [1])
      panic_exit ("a/b mismatch for label %d", index);
    int even_odd = egauss [index].dowker [0] + egauss [index].dowker [1];
    if (even_odd % 2 != 1)
      panic_exit ("label %d not matched to even/odd pair of crossings", index);
    
    // if we get to here, things are probably good
    gauss [egauss [index].dowker [0]] = egauss [index].above_below [0] * egauss [index].dowker [1];
    gauss [egauss [index].dowker [1]] = egauss [index].above_below [1] * egauss [index].dowker [0];
    sign [egauss [index].dowker [0]] = sign [egauss [index].dowker [1]] = egauss [index].sign [0];   
  }
  for (int k = 0; k <= max_label; k++)
    egauss [k].dowker [0] = 0;
}

void get_knot_gauss (char value) {
  static int expect = GAUSS_ab;
  static int above_below;
  static int label;
  static int N = 0;         // local record of number of crossings
  static int max = 0;       // maximum label seen
  int sign2;
  static int linenum = 0;
  static int dowker_label = 1;
  
  switch (value) {
  case 'a':
  case 'b':
    if (expect != GAUSS_ab) 
      panic_exit ("not expecting an 'a' or a 'b'");
    if (value == 'a')
      above_below = GAUSS_ABOVE;
    else
      above_below = GAUSS_BELOW;
    expect = GAUSS_DIGIT;
    label = 0;
    break;
  
  case '|':
  case '/':  // allow this to avoid issues with UNIX pipes
    panic_exit ("reading gauss files for links not yet implemented");
    break;

  case '+':
  case '-':
    if (expect != GAUSS_SIGN)
      panic_exit ("not expecting a '+' or a '-'");

    if (label >= max_xings)   // note: label is allowed to be 0
      panic_exit ("gauss label of %d is out of bounds", label);

    N++;

    if (label > max) 
      max = label;
    
    if (value == '+') 
      sign2 = 1;
    else 
      sign2 = -1;
    add_gauss_crossing (above_below, label, dowker_label++, sign2);

    expect = GAUSS_ab;
    break;

  case '\n':
    linenum++;

    if (expect != GAUSS_ab)
      panic_exit ("new line character in inappropriate location");
    
    if (N % 2 == 1) 
      panic_exit ("found odd number of crossings");

    // NOTE: max is no longer required that max * 2 == N, thereby allowing
    // for gauss codes of the form b108+a102+b333+a108+b102+a333+ which 
    // corresponds to the same knot diagram as does b1+a2+b3+a1+b2+a3+

    // do stuff
    nxings = N / 2;
    extendedGCtoDC (max);
    max = N = 0;
    dowker_label = 1;

    break;

  case '\r':  // ignore stupid carriage returns on Windows
    return;

  default:
    if (!isdigit ((int) value)) 
      panic_exit ("invalid character `%c' in gauss code", value);

    // note: expect not examined here!
    label = 10 * label + (int) (value - '0');
    expect = GAUSS_SIGN;  // might actualy get a digit next
    break;
  }
}

//void invalid_gauss () 
bool get_knot_gauss () {
  char value;
  while (true) {
    if (fscanf (fpin, "%c", &value) == EOF) return false;
    get_knot_gauss (value);
    if (value == '\n') break;
  }

  int n2 = 2 * nxings;

  //printem ();

  for (int i = 1; i <= n2; i++) 
    if (gauss [i] == 0) {
      printem (stderr);
      panic_exit ("not a valid gauss code");
    }

  makeDiagramFromGauss (diag, gauss, sign, nxings);
  return true;
}

bool dont_write_blank_unknots = false;

void put_knot_gauss () {
  if (diag->nxings == 0 && dont_write_blank_unknots)
    fprintf (fpout, "a1-b1-\n");
  else
    GaussCode (fpout, diag);
}

void put_knot_ewing_millett () {
  static bool first_call = true;
  static int num = 1;
  if (first_call) 
    first_call = false;
  else 
    fprintf (fpout, "\n");

  fprintf (fpout, "knot %d\n", num++);

  //if (diag->nxings == 0 && dont_write_blank_unknots)
  if (diag->nxings < 3 && !diag->is_link)
    fprintf (fpout, "1+3c2a2d3d\n2+1b3b3a1c\n3-2c2b1a1d\n");
  else
    EwingMillettCode (fpout, diag);
}  

void put_knot_num_xings () {
  fprintf (fpout, "%d\n", diag->nxings);
}

void put_knot_allgauss () {
  AllGaussCodes (fpout, diag);
}

int N_randomly_chosen_gcs = 0;

void put_knot_Ngauss () {
  NGaussCodes (fpout, diag,  N_randomly_chosen_gcs);
}

void put_knot_jhomfly () {
  if (diag->nxings == 0 && dont_write_blank_unknots)
    fprintf (fpout, "1\n2\n0 1 0 -1\n1\n0 1\n");
  else
    JenkinsHOMFLY (fpout, diag);

}

void put_knot_ignore () {
  // do nothing
}

bool (*get_knot) () = get_knot_three_file;
void (*put_knot) () = put_knot_gauss;

#define usag(X)  if (wiki) fprintf (stderr, " ");fprintf (stderr, X)

void usage (char *error_text) {
  if (strlen (error_text) > 4)
    fprintf (stderr, "*** error: %s\n\n", error_text);

  bool wiki = false;
  if (!strcmp (error_text, "wiki")) {
    wiki = true;
    fprintf (stderr, " ");
  }


  fprintf (stderr, "xinger version %s %s\n \n", __DATE__, __TIME__);
  usag ("Usage: xinger [options/operators] input-xing-sequence output-xing-sequence\n");
  usag ("       or\n");
  usag ("       xinger [options/operators] input-xing-sequence > output-xing-sequence\n");
  usag ("       or\n");
  usag ("       some-xing-producing-program | xinger -- [options/operators] > output-xing-sequence\n");
  usag ("       Program reads a sequence of knots with extended gauss or dowker codes (dowker codes\n");
  usag ("       plus orientation information),\n");
  usag ("       then performs a (possibly null) series of operations on the crossings in the knot,\n");
  usag ("       then writes out \n");
  usag ("       \n");
  usag ("       Can read / write:\n");
  usag ("                   --- old style triple file format (wastes disc space)\n");
  usag ("                   --- compact binary format\n");
  usag ("                   --- gauss code format (default format for writing)\n");
  usag ("       \n");
  usag ("       \n");
  usag ("       xing operators:\n");
  usag ("       +r          remove crossings using R1 and R2 moves until minimized\n");
  usag ("       +R          remove crossings using R1 and R2 moves until minimized, ignoring LIMIT\n");
  usag ("       +1          remove crossings using only R1\n");
  usag ("       +x          pick a crossing uniformly at random and flip it\n");
  usag ("       +c          pick a clasp uniformly at random and flip it\n");
  usag ("       +p          pick a positive crossing uniformly at random and flip it\n");
  usag ("       +n          pick a negative crossing uniformly at random and flip it\n");
  usag ("       +P          pick a positive super crossing uniformly at random and flip it\n");
  usag ("       +N          pick a negative super crossing uniformly at random and flip it\n");
  usag ("       +i          pick the first isolated clasp (if it exists) and flip it\n");
  usag ("       +o          reverse the orientation of the knot\n");
  usag ("       +z          reflect the knot in the z-direction\n");
  usag ("       +Z          flip all crossings randomly\n");
  usag ("       +b          pick a crossing uniformly at random and recombo it\n");
  usag ("       +S          pick a super-crossing uniformly at random and recombo it\n");
  usag ("       +l          find standard labellings for crossings\n");
  usag ("       +k          write number of crossings (useful together with the -wn option)\n");
  usag ("       +K          write number of isolated clasps (useful together with the -wn option)\n");
  usag ("       +D          write debug / diagnostic information about each diagram\n");
  usag ("       -prog PROG  run program PROG (not implemented yet)\n");
  usag ("       \n");
  usag ("       \n");
  usag ("       xing options:\n");
  usag ("       -p PROB     set probability value to PROB (default is 1)\n");
  usag ("       -L LIMIT    minimum number of crossings during removal (default is 0, ignored using -R)\n");
  usag ("        \n");
  usag ("       \n");
  usag ("       read / write options:\n");
  usag ("       -r3         read old style from three files\n");
  usag ("       -rg         read gauss format (only valid with the `--' option, ignored otherwise)\n");
  usag ("       -rb         read binary format only valid with the `--' option, ignored otherwise)\n");
  usag ("       -rx         read KnotPlot X-file format (only valid with the `--' option, ignored otherwise)\n");
  usag ("       -wg         write a gauss code (the default)\n");
  usag ("       -wag        write all gauss codes consistent with knot shadow, ignoring chirality\n");
  usag ("       -wag2       write all gauss codes consistent with knot shadow, considering chirality\n");
  usag ("       -wNg N      write N gauss codes, choosing signs of crossing randomly\n");
  usag ("       -w3         write old style to three files, wasting disc space\n");
  usag ("       -wb         write binary\n");
  usag ("       -wj         write Jenkins HOMFLY format (default style)\n");
  usag ("       -wjo        write Jenkins HOMFLY format (old style)\n");
  usag ("       -wjn        write Jenkins HOMFLY format (new style)\n");
  usag ("       -wkg        write knot group representation\n");
  usag ("       -wem        write Ewing / Millett format\n");
  usag ("       -wn         write nothing (no knots), useful for checking for errors\n");
  usag ("       -wx         write number of crossings\n");
  usag ("       -anc        write average number of isolated clasps\n");
  usag ("       -nx MIN MAX only write knots that have number crossing >= MIN and <= MAX\n");
  usag ("       -v          produce verbose output\n");
  usag ("       -q          run quietly (the default)\n");
  usag ("       -s SEED     set random number seed to SEED (use this only for testing purposes)\n");
  usag ("       -f FORMAT   set 3-file format base name to FORMAT\n");
  usag ("       -zok        OK to output knots with zero crossings\n");
  usag ("       -znb        don't write blank lines for knots with zero crossings\n");
  usag ("       --          read data from standard input instead of opening a file\n");
  usag ("       \n");
  usag ("       for complete information see Wiki page http://ewok.sfsu.edu/wiki/index.php/Xinger\n");

  exit (1);
}

#define MAX_PROGRAM_LENGTH 24


void add_operator (char *prog, char com) {
  if (com == (char) 0)
    panic_exit ("missing operator or unknown option");
  int len = strlen (prog);
  if (len == MAX_PROGRAM_LENGTH)
    panic_exit ("program too long (maximum length is %d)", MAX_PROGRAM_LENGTH);
  prog [len] = com;
}

char operators [NUM_COMMANDS] = {  // not used yet
  REIDREM_via_R1R2,
  REIDREM_via_R1,
  FLIP_XING,
  CLASP_SP,
  PRINT_ISO_CLASPS,
  REVERSE_ORIENTATION,
  REFLECT_IN_Z};


int main (int argc, char **argv) {
  bool verbose = false;
  bool writing_nothing = false;

  if (argc < 2) 
    usage ("");

  char *input_file = (char *) NULL;
  char *output_file = (char *) NULL;

  int index = 1;
  int RNGoffset = 0, RNGseed;
  time_t the_time;
  time (&the_time); 
  RNGseed = (int) the_time;
  char *program = (char *) calloc (MAX_PROGRAM_LENGTH + 2, sizeof (char));
  double prob = 1.0;
  bool zok = false;
  bool read_from_stdin = false;
  bool write_KnotPlot_script_and_scenery = false;
  bool write_anc = false;
  int min_nxings = 0;   // not to be confused with the variable below
  int min_nxings_write = 0;
  int max_nxings_write = HUGE_NUMBER;
  //fprintf (stderr, "xinger: test version\n"); fflush (stderr);
  bool output_junk = false;
  
  while (index < argc && (argv [index][0] == '-' || argv [index][0] == '+')) {
    if (!strcmp (argv [index], "-r3")) {
      get_knot = get_knot_three_file;
      input_format = THREE_FILE;
    }
    else if (!strcmp (argv [index], "-w3")) {
      put_knot = put_knot_three_file;
      output_format = THREE_FILE;
    }
    else if (!strcmp (argv [index], "-wg")) {
      put_knot = put_knot_gauss;
      output_format = GAUSS;
    }
    else if (!strcmp (argv [index], "-wag") || !strcmp (argv [index], "-wag2")) {
      put_knot = put_knot_allgauss;
      output_format = ALLGAUSS;
      if (strlen (argv [index]) == 5) 
	  wag2 = true;
    }
    else if (!strcmp (argv [index], "-wNg")) {
      if (++index >= argc) 
	panic_exit ("expecting an argument for -wNg");
      if (!alldigits (argv [index]))
	panic_exit ("expecting an integer argument for -wNg");
      put_knot = put_knot_Ngauss;
      output_format = NGAUSS;
      N_randomly_chosen_gcs = atoi (argv [index]);
    }
    else if (!strcmp (argv [index], "-rb")) {
      get_knot = get_knot_binary;
      input_format = BINARY;
    }
    else if (!strcmp (argv [index], "-rg")) {
      get_knot = get_knot_gauss;
      input_format = GAUSS;
    }
    else if (!strcmp (argv [index], "-rx")) {
      get_knot = get_knot_KnotPlot_X;
      input_format = KNOTPLOT_X;
    }
    else if (!strcmp (argv [index], "--")) {
      read_from_stdin = true;
    }
    else if (!strcmp (argv [index], "-wb")) {
      put_knot = put_knot_binary;
      output_format = BINARY;
    }
    else if (!strcmp (argv [index], "-wem")) {
      put_knot = put_knot_ewing_millett;
      output_format = EMCODE;
    }
    else if (!strcmp (argv [index], "-wj")) {
      put_knot = put_knot_jhomfly;
      output_format = JHOMFLY;
      //Jenkins_old_style = false;
    }
    else if (!strcmp (argv [index], "-wjo")) {
      put_knot = put_knot_jhomfly;
      output_format = JHOMFLY;
      Jenkins_old_style = true;
    }
    else if (!strcmp (argv [index], "-wjn")) {
      put_knot = put_knot_jhomfly;
      output_format = JHOMFLY;
      Jenkins_old_style = false;
    }
    else if (!strcmp (argv [index], "-wkg")) {
      put_knot = put_knot_kgroup;
      output_format = KGROUP;
    }
    else if (!strcmp (argv [index], "-wn")) {
      put_knot = put_knot_ignore;
      output_format = IGNORE;
      writing_nothing = true;
    }
    else if (!strcmp (argv [index], "-wx")) {
      put_knot = put_knot_num_xings;
      output_format = NUM_XINGS;
    }
    else if (!strcmp (argv [index], "-colour")) {
      colour_it = true;
      output_junk = true;
    }
    else if (!strcmp (argv [index], "-v")) {
      verbose = true;
    }
    else if (!strcmp (argv [index], "-f")) {
      panic_exit ("option `-f' not implemented yet");
    }
    else if (!strcmp (argv [index], "-p")) {
      if (index == argc - 1)
        usage ("expecting a value for the clasp strand passage probability");
      prob = atof (argv [++index]);
    }
    else if (!strcmp (argv [index], "-L")) {
      if (index == argc - 1)
        usage ("expecting a value for the minimum number of crossings");
      min_nxings = atoi (argv [++index]);
    }
    else if (!strcmp (argv [index], "-q")) {
      verbose = false;
    }
    else if (!strcmp (argv [index], "-s")) {
      if (index == argc - 1)
        usage ("expecting a value for random number seed");
      RNGseed =  atoi (argv [++index]);
      fprintf (stderr, " *** WARNING: setting random number seed to %d for testing purposes!\n", RNGseed);
    }
    else if (!strcmp (argv [index], "-so")) {
      if (index == argc - 1)
        usage ("expecting a value for random number seed offset");
      RNGoffset =  atoi (argv [++index]);
    }
    else if (!strcmp (argv [index], "-zok")) {
      zok =  true;
    }
    else if (!strcmp (argv [index], "-anc")) {
      write_anc =  true;
    }
    else if (!strcmp (argv [index], "-znb")) {
      dont_write_blank_unknots =  true;
    }
    else if (!strcmp (argv [index], "-kps")) {
      write_KnotPlot_script_and_scenery = show_iso_clasps = true;
      output_format = KNOTPLOT_S;
    }
    else if (!strcmp (argv [index], "-xkps")) {
      write_KnotPlot_script_and_scenery = true;
      show_iso_clasps = false;
      output_format = KNOTPLOT_S;
    }
    else if (!strcmp (argv [index], "-ray")) {
      fp_ray = open_file ("rays.dat", "w");
      //index++;
    }
    else if (!strcmp (argv [index], "-wiki")) {
      usage ("wiki");
    }
    else if (!strcmp (argv [index], "-nx")) {
      if (index + 2 >= argc) 
	panic_exit ("expecting MIN and MAX values");
      min_nxings_write = atoi (argv [++index]);
      max_nxings_write = atoi (argv [++index]);
    }
    else {
      if (argv [index][0] == '-')
	warn_deprecated = true;
      add_operator (program, argv [index][1]);
    }
    index++;
  }

  //  add_operator (program, 'r'); fprintf (stderr, " ***** WARNING, temporary fix, tell Rob to remove this!\n\n");

  if (input_format == THREE_FILE && read_from_stdin) 
    panic_exit ("can't use three-file format when reading from standard input, use either `-rb' or `-rg'");

  if (input_format == UNKNOWN && read_from_stdin) {
    get_knot = get_knot_gauss;
    input_format = GAUSS;
  }

  if (read_from_stdin) {
    if (input_format == KNOTPLOT_X) 
      panic_exit ("option `-rx' together with option `--' not yet implemented");

    if (argc != index)
      usage ("too many arguments on command line");
    allocate_string (&input_file, "stdin");
  }
  else {
    if (!(index == argc - 2 || index == argc - 1)) 
      usage ("inappropriate number of file names given");
    allocate_string (&input_file, argv [index++]);
  }

  if (warn_deprecated) 
    fprintf (stderr, "\n *** specifying xing operators using `-' is deprecated.  use `+' instead\n");

  sRandSimple (RNGseed + RNGoffset);

  char *dfname, *cfname, *ofname;
  char *odfname, *ocfname, *oofname;

  switch (input_format) {
  case THREE_FILE:
    cfname = (char *) calloc (strlen (input_file) + 16, sizeof (char));
    dfname = (char *) calloc (strlen (input_file) + 16, sizeof (char));
    ofname = (char *) calloc (strlen (input_file) + 16, sizeof (char));
    sprintf (cfname, "crossings_%s.txt", input_file);
    sprintf (dfname, "dowker_%s.txt", input_file);
    sprintf (ofname, "orient_%s.txt", input_file);
    fpic = open_file (cfname, "r");
    fpid = open_file (dfname, "r");
    fpio = open_file (ofname, "r");
    if (verbose && index != argc || writing_nothing)
      printf ("reading from %s, %s and %s\n", cfname, dfname, ofname);
    break;
  default:
    int magic = 0;
    unsigned char version;

    if (read_from_stdin) {
      fpin = stdin;
      switch (input_format) {
      case BINARY:
#ifdef WIN32
        _setmode (_fileno (stdin), _O_BINARY);
#endif
        read_int (&magic, fpin);
        if (magic != MAGIC_XING)
          panic_exit ("expecting a XING file");
        get_knot = get_knot_binary;
        fread (&version, 1, 1, fpin);
        if (version != 0) 
          panic_exit ("unknown XING file version");
        break;
      case GAUSS:
        get_knot = get_knot_gauss;
        egauss = (ExtendedGaussCode *) calloc (max_xings + 1, sizeof (ExtendedGaussCode));
        break;
      default:
        panic_exit ("the impossible has happened");
      }
    }
    else {
      fpin = open_file (input_file, "rb");

      // try to understand what we're dealing with here
      read_int (&magic, fpin);
      if (magic == MAGIC_XING) {
        input_format = BINARY;
        get_knot = get_knot_binary;
        fread (&version, 1, 1, fpin);
        if (version != 0) 
          panic_exit ("unknown XING file version");
      }
      else {
        char ch = (magic >> 24) & 0xff;
        //printf ("ch is %c\n", ch);
        if (ch == 'a' || ch == 'b') {
          input_format = GAUSS;
          get_knot = get_knot_gauss;
          egauss = (ExtendedGaussCode *) calloc (max_xings + 1, sizeof (ExtendedGaussCode));
        }
        else if (isdigit ((int) ch)) {
          input_format = KNOTPLOT_X;
          get_knot = get_knot_KnotPlot_X;
          coords = (Coord *) calloc (2 * max_xings + 1, sizeof (Coord));
        }
        else 
          panic_exit ("unknown input file format");
   
        
        fseek (fpin, 0, SEEK_SET);  // rewind file to start
      }
    }
    break;
  }

  FILE *fpinfo;

  if (index == argc) {
    if (output_format == THREE_FILE)
      usage ("can't output three files to standard output!");
    if (output_format == KNOTPLOT_S)
      usage ("can't output KnotPlot script and scenery to standard output!");
    fpout = stdout;
    fpinfo = stderr;
#ifdef WIN32
    _setmode (_fileno (stdout), _O_BINARY);
#endif
    allocate_string (&output_file, "stdout");
  }
  else {
    allocate_string (&output_file, argv [index++]);
    switch (output_format) {
    case THREE_FILE:
      if (input_format == THREE_FILE && !strcmp (input_file, output_file))
        panic_exit ("input and output file names are the same!");

      ocfname = (char *) calloc (strlen (output_file) + 16, sizeof (char));
      odfname = (char *) calloc (strlen (output_file) + 16, sizeof (char));
      oofname = (char *) calloc (strlen (output_file) + 16, sizeof (char));
      sprintf (ocfname, "crossings_%s.txt", output_file);
      sprintf (odfname, "dowker_%s.txt", output_file);
      sprintf (oofname, "orient_%s.txt", output_file);
      fpoc = open_file (ocfname, "w");
      fpod = open_file (odfname, "w");
      fpoo = open_file (oofname, "w");
      if (verbose)
        printf ("writing to %s, %s and %s\n", ocfname, odfname, oofname);
      break;
    default:  
      fpout = open_file (output_file, "wb");   // NOTE: need to specify binary (b) mode in order
                                               // to work properly on Windows 
      if (verbose)
        printf ("writing to %s\n", output_file);
    }
  }

  if (fpout != stdout && verbose) { 
    printf ("running program `%s' on the knots\n", program);
  }

  diag = initDiagram (max_xings, prob, zok);
  diag->min_nxings = min_nxings;
  diag->max_nxings_write = max_nxings_write;
  diag->min_nxings_write = min_nxings_write;

  gauss = (int *) calloc (2 * max_xings + 1, sizeof (int));
  sign = (int *) calloc (2 * max_xings + 1, sizeof (int));

  if (write_KnotPlot_script_and_scenery) {
    if (get_knot ()) 
      output_KnotPlot_script_and_scenery (output_file, diag);
    exit (102);
  }

  int total_anc = 0;

  while (get_knot ()) {
    executeProgram (diag, program, fpout);
    if (diag->nxings >= diag->min_nxings_write &&
	diag->nxings <= diag->max_nxings_write) {
      put_knot ();
      if (write_anc) 
	total_anc += findIsolatedClasps (diag);
    }
    knotnumber++;
  }

  if (colour_it) 
    output_colours (diag);

  if (write_anc) 
    printf ("%f\n", (double) total_anc / (double) (knotnumber - 1));

  switch (output_format) {
  case THREE_FILE:
    fprintf (fpoc, "\n");
    break;
  case KGROUP:
    fprintf (fpout, "END FILE\n");
    break;
  }

  if (fpout != stdout && verbose) {
    printf ("%d isolated clasps\n%d isolated clasp flips\n%d normal flips\n%d crossing removals\n", 
      diag->iso_clasps, diag->iso_clasp_flips, diag->xing_flips, diag->removals);
    if (diag->min_nxings != 0) 
      printf ("limited xing removal to %d\n", diag->min_nxings);
  }

  if (fp_ray)
    fclose (fp_ray);

  return 0;
}

