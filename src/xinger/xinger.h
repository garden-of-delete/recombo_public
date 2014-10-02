#ifndef XINGER_H
#define XINGER_H

typedef struct {
  double x;
  double y;
  double z;
  double angle;
} Coord;

#define RAD_TO_DEG   57.2957795130823208768    /*  (180 / PI)  */
#define DEG_TO_RAD   0.0174532925199432957692  /*  (PI / 180)  */


// the following are powers of 2 (i.e., single bit values) to allow for the possibility
// that the software will (incorrectly) classify a xing as more than one type
// see the function classifyXings() 

#define XING_KIND_NORMAL        0
#define XING_KIND_NUGATORY      1
#define XING_KIND_ISOCLASP      2
#define XING_KIND_SUPERXING     4

typedef struct _xing {
  struct _xing *xtop [2];
  struct _xing *xbot [2];
  struct _xing *xother;    // can be used for various things
  int stop [2];
  int sbot [2];
  int lab [2];      // lab [0] is crossing label associated with strand 0 
  int sign;
  int ID;
  int mark;         // used to mark stuff
  int location;     // location in diag->xing [] array
  int start;  // strand to start with later, either 0 or 1
  int gauss;
  int last_strand;
  double x, y, z, angle;
  int test;
  int region [4];
  int kind;
  char emtop [2];   // Ewing / Millett labels for top of strands
  char embot [2];
  char *colour;
} Xing;

typedef Xing *XingPtr;

typedef struct {
  int nxings;
  int max_xings;
  int *loc;
  int *scratch;
  int min_nxings;
  XingPtr *xing;
  XingPtr first_xing;
  XingPtr last_xing;
  double p;
  int xing_flips;
  int iso_clasps;
  int iso_clasp_flips;
  int num_knots_with_iso_clasps;
  int removals;
  bool zok;
  bool is_link;
  int nregions;
  char *gauss_code;
  char *cscratch;
  int min_nxings_write;
  int max_nxings_write;
} Diagram;

typedef Diagram *DiagramPtr;

DiagramPtr initDiagram (int max_xings);
DiagramPtr initDiagram (int max_xings, double p, bool non_zero_nxings);


void makeDiagramFromGauss (DiagramPtr diag, int *gauss, int *sign, int nxings);
void makeDiagramFromGauss (DiagramPtr diag, int *gauss, int *sign, Coord *coords, int nxings);

void makeGaussFromDiagram (int *gauss, int *sign, int nxings, DiagramPtr diag);
  
bool removeXingsR1 (DiagramPtr diag);
bool removeXingsR2 (DiagramPtr diag);
int findIsolatedClasps (DiagramPtr diag);
int printIsolatedClasps (DiagramPtr diag);
bool flipIsolatedClasp (DiagramPtr diag);

void executeProgram (DiagramPtr diag, char *program, FILE *fp);

#define NUM_COMMANDS 7

#define REIDREM_via_R1R2    'r'
#define REIDREM_via_R1R2_0  'R'
#define REIDREM_via_R1      '1'
#define FLIP_XING           'x'
#define CLASP_SP            'c'
#define FLIP_POS_XING       'p'
#define FLIP_NEG_XING       'n'
#define FLIP_POS_SUPERXING  'P'
#define FLIP_NEG_SUPERXING  'N'
#define PRINT_ISO_CLASPS    'I'
#define REVERSE_ORIENTATION 'o'
#define REFLECT_IN_Z        'z'
#define CHECK_DIAGRAM       'C'
#define RECOMBO             'b'
#define SUPER_FLIP          'S'
#define WRITE_NUM_XINGS     'k'
#define WRITE_NUM_CLASPS    'K'
#define STANDARD_LABELLING  'l'
#define FLIP_ALL_XINGS      'Z'
#define DETECT_R3_R2_SIMPS  '3'
#define PRINT_ENTANGLING    'e'
#define PRINT_BRIDGE_INFO   'B'
#define PRINT_DEBUG_INFO    'D'
#define FLIP_FIRST_ISOCLASP 'i'

#define MAGIC_XING 1481199175  // word "XING" as an big-endian int
void DTcode (FILE *fpb, DiagramPtr diag);                        // for binary output
void DTcode (FILE *fpc, FILE *fpd, FILE *fpo, DiagramPtr diag);  // for plain text output
void GaussCode (FILE *fp, DiagramPtr diag);
void AllGaussCodes (FILE *fp, DiagramPtr diag);
void NGaussCodes (FILE *fp, DiagramPtr diag, int N);
void JenkinsHOMFLY (FILE *fp, DiagramPtr diag);
void KnotGroup (FILE *fp, DiagramPtr diag);
void EwingMillettCode (FILE *fp, DiagramPtr diag);

void output_KnotPlot_script_and_scenery (char *basename, DiagramPtr diag);

void standardLabelling (DiagramPtr);
void findStandardLabelling (DiagramPtr);
void checkDiagram (DiagramPtr, char *);
bool labelDiagram (DiagramPtr);
bool labelDiagram2 (DiagramPtr);
//bool labelDiagram3 (DiagramPtr, bool);
bool isNugatory (XingPtr xp);

extern FILE *fp_ray;
extern bool superv;
extern bool colour_it;
extern bool show_iso_clasps;
extern bool wag2;
extern bool Jenkins_old_style;

void output_colours (DiagramPtr diag);

void KnotGroupRepresentation (DiagramPtr diag);
#define PIBY2        1.57079632679489661923
#define HUGE_NUMBER  2000000000


#endif
