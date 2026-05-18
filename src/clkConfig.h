#pragma once

#include "clkTypes.h"

void bfacf_set_probabilities(ComponentCLKPtr comp, double pm2, double p0, double pp2);
void bfacf_set_probabilities(ComponentCLKPtr comp, double z);
void bfacf_set_probabilities(CubicLatticeKnotPtr knot, double z);
void bfacf_lattice_info(CubicLatticeKnotPtr knot);
void pexit(char* s);
