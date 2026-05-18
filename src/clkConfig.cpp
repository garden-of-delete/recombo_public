#include "clkConfig.h"
#include "clkConstants.h"

#include <cstdio>
#include <cstdlib>

void bfacf_set_probabilities(ComponentCLKPtr comp, double pm2, double p0, double pp2)
{
   comp->p_minus2 = pm2;
   comp->p_0 = p0;
   comp->p_plus2 = pp2;

   // probably doesn't make much of a difference,
   // but precompute the following anyway:

   comp->p_4p2 = 4.0 * comp->p_plus2;
   comp->p_03p2 = comp->p_0 + 3.0 * comp->p_plus2;
   comp->p_m23p2 = comp->p_minus2 + 3.0 * comp->p_plus2;
   comp->p_2p0 = 2.0 * comp->p_0;
   comp->p_2p02p2 = 2.0 * comp->p_0 + 2.0 * comp->p_plus2;
}

void bfacf_set_probabilities(ComponentCLKPtr comp, double z)
{
   double p_plus2 = (z * z) / (1.0 + 3.0 * z * z);
   double p_0 = (1.0 + z * z) / (2.0 * (1.0 + 3.0 * z * z));
   double p_minus2 = 1.0 / (1.0 + 3.0 * z * z);
   comp->z = z;
   bfacf_set_probabilities(comp, p_minus2, p_0, p_plus2);
}

void bfacf_set_probabilities(CubicLatticeKnotPtr knot, double z)
{
   ComponentCLKPtr comp = knot->fcomp;
   while (comp)
   {
      bfacf_set_probabilities(comp, z);
      comp = comp->next;
   }
}

void bfacf_lattice_info(CubicLatticeKnotPtr knot)
{
   int min, max;
   int minloc, maxloc;
   minloc = maxloc = 0;
   min = max = knot->alt_lattice [0];
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (knot->alt_lattice [i] > max)
      {
         max = knot->alt_lattice [i];
         maxloc = i;
      }
      if (knot->alt_lattice [i] < min)
      {
         min = knot->alt_lattice [i];
         minloc = i;
      }
   }
   printf("%d, min value of %d at %d, max value of %d at %d\n", knot->nedges_total, min, minloc, max, maxloc);
}

void pexit(char *s)
{
   fprintf(stderr, "%s\n", s);
   fflush(stderr);
   exit(102);
}
