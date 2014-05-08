#ifndef LEGACYBFACF_H
#define LEGACYBFACF_H

// from bfacf.h
#ifdef _WIN32   // disable annoying warnings
#pragma warning (disable: 4996)
#endif

#include "legacy.h"

typedef struct {
  double p_plus2;
  double p_0;
  double p_minus2;
  double p_4p2;
  double p_03p2;
  double p_m23p2;
  double p_2p0;
  double p_2p02p2;
} qProb;
 
typedef struct _edge {
  struct _edge *prev;   // pointer to previous edge in doubly linked list of edges
  struct _edge *next;   // pointer to next edge in doubly linked list of edges
  int dir;              // direction of edge, integer from 0 to 5, inclusive
  int increment [3];    // increment corresponding to dir above 
  int start [3];        // position of start of edge, a 3-tuple of integers
  int locpool;          // location of this edge in edge pool
  bool frozen;          // whether is frozen
  void  *comp;          // component this edge is located in
  int colour;           // edge colour
  int ID;               // unique ID (for debugging purposes), assigned when edge created
  int scratch;          // temporary space 
} Edge;

typedef Edge *EdgePtr;

typedef struct {
  EdgePtr ep1, ep2;
} EdgePair;

typedef EdgePair *EdgePairPtr;

#define COMPONENT_CLK_FLAG_OPEN  1  // open-ended

typedef struct _comp {
  int ID;
  int flags;             // properties of this component
  int nedges;            // number of edges in knot component
  int minedges;          // minimum number of edges permitted in this component
  int maxedges;          // maximum number of edges permitted in this component
  void *clkp;            // which clkp this component is in
  EdgePtr first_edge;    // pointer to first edge in knot component
  EdgePtr last_edge;     // pointer to last edge in knot component
  struct _comp *prev;    // previous component
  struct _comp *next;    // next component

  double z;              // fugacity parameter, might be different for different knots

  // values of the next three parameters computed from z value
  double p_minus2;       // p(-2), probability of doing a -2 move 
  double p_0;            // p(0), probability of doing a 0 move
  double p_plus2;        // p(+2), probability of doing a +2 move

  // the following are precomputed in hope of an improvement in running time
  double p_4p2;          // 4 * p(+2)
  double p_03p2;         // p(0) + 3 * p(+2)
  double p_m23p2;        // p(-2) + 3 * p(+2)  
  double p_2p0;          // 2 * p(0)
  double p_2p02p2;       // 2 * p(0) + 2 * p(+2)
} ComponentCLK;          // cubic lattice knot component

typedef ComponentCLK *ComponentCLKPtr;

typedef struct {
  int loffset [3];       // offset of first vertex in self-avoidance lattice
  int nedges_total;      // total number of edges in knot (note: not number in edge pool)
  int nedges_alt;        // not used except by KnotPlot (for now)

  EdgePtr *edgepool;     // pool of pointers to edges (array of size `poolsize')
  int poolsize;          // maximum number of edges in pool
  int nfrozen;           // number of `frozen edges', frozen edges are never chosen
                         // (this allows for open-ended strings or the inclusion
                         // of topological obstructions)
  bool auto_recentre;    // whether or not to recentre after an edge hit
  bool filter_by_energy;

  char *lattice;         // lattice used for self-avoidance checking
  int  *alt_lattice;     // alternate lattice, allocated on demand

  ivector max_range;     // maximum range allowed

  ComponentCLKPtr fcomp; // first cubic lattice knot component in a linked list
  ComponentCLKPtr lcomp; // last cubic lattice knot component in a linked list
  int ncomps;            // number of components


  // the following are to keep track of various things that may or may not be of interest
  // but are not actually used in the algorithm
  int edge_hits [6];     // number of times knot hits edge of self-avoidance lattice 
  int success_minus2;    // number of successful -2 moves
  int success_0;         // number of successful 0 moves
  int success_plus2;     // number of successful +2 moves
  int p2dir [6];
  int pm2dir [6];
  int total_sampled;
  int total_output;
  int wrong_length;

  FILE *fpout;

  int iteration;         // iteration count

  qProb *q_prob;
  int q_prob_size;

  // related to knot energy
  // double T;              // absolute temperature
  //double kappa;

  // related to Multiple Markov Chain (MMC)
  int swaps, attempted_swaps;

} CubicLatticeKnot;

typedef CubicLatticeKnot *CubicLatticeKnotPtr;



#define Boltzmann_constant  1.3806503e-23
#define Tesi_v             -0.26
#define Z_CRITICAL          0.2134

bool add_edge_to_knot (CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end);
CubicLatticeKnotPtr bfacf_input_start_configuration (FILE *fpin);
CubicLatticeKnotPtr bfacf_input_start_configuration (FILE *fpin, bool ignore_mid);
CubicLatticeKnotPtr bfacf_input_start_configuration (FILE *fp, bool ignore_mid, int poolsize);
void bfacf_init_pool_and_lattice (CubicLatticeKnotPtr clkp, bool ignore_mid, int poolsize, ivector min, ivector max);
bool perform_move (CubicLatticeKnotPtr knot);
bool perform_move_q (CubicLatticeKnotPtr knot);
bool perform_move_kjc (CubicLatticeKnotPtr knot);
bool perform_move_move0_only (CubicLatticeKnotPtr knot);
bool perform_move_pivot (CubicLatticeKnotPtr knot);

int bfacf_save_stick (int nsticks, CubicLatticeKnotPtr clkp);
bool type_minus2_moves_possible (CubicLatticeKnotPtr clkp);
int number_type0_moves (CubicLatticeKnotPtr clkp);
void bfacf_number_moves (CubicLatticeKnotPtr clkp, int &nm2, int &n0, int &np2);
int bfacf_number_sticks (CubicLatticeKnotPtr clkp);
bool freezeEdge (CubicLatticeKnotPtr knot, ivector start);
void delete_Edge (CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep);

void swapEdgesPool (CubicLatticeKnotPtr clkp, int loc1, int loc2);
void freezeEdge (CubicLatticeKnotPtr clkp, EdgePtr ep);
void freezeEdge (EdgePtr ep);
void thawEdge (CubicLatticeKnotPtr clkp, EdgePtr ep);
void thawEdge (EdgePtr ep);

extern int clk_minedges, clk_maxedges;
void clk_set_nedge_limit (int min, int max);
void clk_set_nedge_limit (CubicLatticeKnotPtr clkp, int min, int max, int ID);
EdgePtr clk_get_edge (CubicLatticeKnotPtr clkp, ivector start);

#ifdef CLAMP
#undef CLAMP
#endif
#define CLAMP(V,L,H)   V = MAX (L, MIN (V, H))

extern int opposite [6];
extern ivector increment_NEWSUD [6];
extern int turn [6][4];
extern bool perp [6][6];
extern bool anti [6][6];
extern char clk_dir_name [6][6];
extern int kross [6][6];

#define MOVE_NORTH 0
#define MOVE_EAST  1
#define MOVE_WEST  2
#define MOVE_SOUTH 3
#define MOVE_UP    4
#define MOVE_DOWN  5

#define MOVE_INVALID -1

#define TWO_BILLION 2000000000
#define HUGE_NUMBER TWO_BILLION

#define DEFAULT_z             0.20815       
#define DEFAULT_WARMUP        1000000
#define DEFAULT_CHECK         5000
#define DEFAULT_ITERATIONS    60000000
#define DEFAULT_SAVE_INTERVAL 5000
#define DEFAULT_POOLSIZE      108000
//#define DEFAULT_POOLSIZE      5080000

#define RED     1
#define GREEN   2
#define YELLOW  (RED|GREEN)   // 3
#define BLUE    4
#define MAGENTA (RED|BLUE)    // 5
#define CYAN    (GREEN|BLUE)  // 6

// from clk_util.h

int bfacf_perform_recombination (CubicLatticeKnotPtr clkp, void mark (EdgePtr, EdgePtr));
int bfacf_perform_recombination (CubicLatticeKnotPtr clkp);
bool perform_recombination(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2);
void bfacf_set_probabilities (ComponentCLKPtr comp, double pm2, double p0, double pp2);
void bfacf_set_probabilities (CubicLatticeKnotPtr knot, double z);
void bfacf_set_probabilities(ComponentCLKPtr comp, double z);
void clk_get_extent (ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot);
bool clk_check_increment (ivector incr);
bool clk_check_lattice (CubicLatticeKnotPtr knot);
void clk_check_increments (CubicLatticeKnotPtr knot);
void clk_par ();
bool clk_check_for_edge_hits (CubicLatticeKnotPtr knot, int dir, ivector test_location);
int clk_direction (ivector incr);
int clk_direction (ivector end, ivector start);
void clk_get_statistics (CubicLatticeKnotPtr clkp, int &n1, int &n2, int &n3, int &n4);
void get_NEWS (char *news, CubicLatticeKnotPtr clkp);
void init_lattice (CubicLatticeKnotPtr clkp, ivector mid);
void fill_lattice (CubicLatticeKnotPtr clkp, ivector min, ivector max);
void fill_lattice (CubicLatticeKnotPtr clkp);
void clear_lattice (CubicLatticeKnotPtr clkp);
void set_lattice (CubicLatticeKnotPtr clkp, char value);
void invert_lattice (CubicLatticeKnotPtr clkp);
void dilate_lattice (CubicLatticeKnotPtr clkp);
void clear_lattice (CubicLatticeKnotPtr clkp, ivector min, ivector max);
bool recentre_knot_in_lattice (CubicLatticeKnotPtr knot);
bool findz (int target_length, char *knotname, double &z);
void copy_array_to_CLKP (ivector *coord, int nedges, CubicLatticeKnotPtr clkp);
void copy_CLKP_to_array (ivector *coord, int &nedges, CubicLatticeKnotPtr clkp);
bool bfacf_strand_pass (CubicLatticeKnotPtr clkp); 
bool bfacf_parsite_pass (CubicLatticeKnotPtr clkp); 
void clk_validate (CubicLatticeKnotPtr clkp, char *s); 
int clk_validate(CubicLatticeKnotPtr clkp);
bool clk_verbose_knot_path_info (CubicLatticeKnotPtr clkp, char *filename);
extern bool clk_recombo_limit, clk_recombo_direct;
extern int clk_min_arc_length_distance;

#define LATTICE_SIZE        256  
#define LATTICE_SLICE_SIZE  (LATTICE_SIZE * LATTICE_SIZE)
#define LATTICE_TOTAL_SIZE  (LATTICE_SIZE * LATTICE_SIZE * LATTICE_SIZE)

#define lat(P)   (P [2] * LATTICE_SLICE_SIZE + P [1] * LATTICE_SIZE + P [0])

#define NEXTOCCUPIED 2   // lattice cell is next to an occupied cell
#define OCCUPIED     1
#define EMPTY        0

#define ARC_LENGTH_DISTANCE_INFINITE 1000000000
#define ARC_LENGTH_DISTANCE_NO_EDGE1 -1
#define ARC_LENGTH_DISTANCE_NO_EDGE2 -2

int clk_arc_length_distance (CubicLatticeKnotPtr clkp, ivector s1, ivector s2);
int clk_arc_length_distance (EdgePtr ep1, EdgePtr ep2);
void energy_report (CubicLatticeKnotPtr, char *);

double energy (ivector *currentKnot, int length);
double energy (ivector A, ivector B);

extern double T;

#define B  1.3806503e-23    // the Boltzmann constant in units of J K^-1

void energy_set_nu (double value);
void energy_set_A (double value);
void energy_set_kappa (double value);
void energy_set_T (double value);

void set_recombo_direct (bool);
double energy (ivector location, CubicLatticeKnotPtr clkp, EdgePtr beg, EdgePtr end);

#define DEFAULT_KAPPA 0.32
#define DEFAULT_KELVINS 310.0
#define DEFAULT_A 0.01
#define DEFAULT_NU  -0.26

void energy_blurt_values (char *);

void save_knot_energy (FILE *fp, CubicLatticeKnotPtr clkp);
void clk_bounding_box (int &X, int &Y, int &Z, CubicLatticeKnotPtr knot);
void bfacf_initialize_expd ();
void MMC_prob_report (char *knot);
bool MMC_swap (double &prob, double zi, int lengthi, double zip1, int lengthip1);
char *get_known_lengths (char *knot);

void  bfacf_set_lattice_sphere (CubicLatticeKnotPtr clkp, int what, double radius, double xcen, double ycen, double zcen, char *blurt);
int bfacf_info_filled (CubicLatticeKnotPtr knot);
bool clkp_always_turns (CubicLatticeKnotPtr clkp);

int numb_tight_clasps (CubicLatticeKnotPtr clkp);
int numb_tight_clasps (CubicLatticeKnotPtr clkp, void mark (ivector, ivector));
bool has_tight_clasp (CubicLatticeKnotPtr clkp);
EdgePtr is_tight_clasp (CubicLatticeKnotPtr clkp, EdgePtr ep);

int numb_parsites (CubicLatticeKnotPtr clkp);
int numb_parsites (CubicLatticeKnotPtr clkp, void mark (EdgePtr, EdgePtr));
bool has_parsite (CubicLatticeKnotPtr clkp);
EdgePtr is_parsite (CubicLatticeKnotPtr clkp, EdgePtr ep);

bool clk_allocation_alt_lattice (CubicLatticeKnotPtr clkp);
void init_probabilities_q (CubicLatticeKnotPtr clkp, int ni, double z, double q, bool blurt);

void print_zknown (bool);
void clk_fill_path_alt (CubicLatticeKnotPtr clkp, bool fill) ;
//extern qProb *q_prob;
//extern int q_prob_size;

// from findz.h


#endif