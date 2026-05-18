// Core types and constants for the BFACF algorithm
#include "clkTypes.h"
#include "clkConstants.h"
#include "clkTables.h"

// Function declarations from refactored modules
#include "clkInit.h"
#include "clkConfig.h"
#include "clkGeometry.h"
#include "clkValidation.h"
#include "clkLattice.h"
#include "clkTopology.h"
#include "bfacfMove.h"
#include "edgePool.h"
#include "clkRecombination.h"

// from bfacf.h
#ifdef _WIN32   // disable annoying warnings
#pragma warning (disable: 4996)
#endif
