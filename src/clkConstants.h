#pragma once

#define Z_CRITICAL          0.2134

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

#define RED     1
#define GREEN   2
#define YELLOW  (RED|GREEN)   // 3
#define BLUE    4
#define MAGENTA (RED|BLUE)    // 5
#define CYAN    (GREEN|BLUE)  // 6

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

#define CLAMP(V,L,H)   V = MAX (L, MIN (V, H))
