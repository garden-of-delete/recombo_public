#include "clkValidation.h"
#include "clkConstants.h"
#include "clkTables.h"
#include "clkGeometry.h"

#include <cstdio>
#include <cstdlib>

int clk_validate(CubicLatticeKnotPtr clkp)
{
    for (int i = 0; i < clkp->nedges_total; i++)
        if (clkp->edgepool[i]->locpool != i)
            return 1;

    int kount = 0;
    ComponentCLKPtr comp = clkp->fcomp;
    while (comp)
    {
        EdgePtr start;
        EdgePtr ep = start = comp->first_edge;
        if (comp->nedges < 4) return 7;
        ivector loc;
        copy_ivector(loc, ep->start);
        for (int i = 0; i < comp->nedges; i++)
        {
            if (ep->prev->next != ep)
                return 2;
            if (ep->next->prev != ep)
                return 3;
            if ((ComponentCLKPtr)ep->comp != comp)
                return 5;
            if (clkp->lattice[lat(ep->start)] != OCCUPIED)
                return 8;
            add_ivector(loc, loc, ep->increment);
            if (!equal_ivector(loc, ep->next->start))
                return 9;
            if (ep->dir != clk_direction(ep->increment))
                return 10;
            ep = ep->next;
        }
        if (ep != start) return 4;
        kount += comp->nedges;
        comp = comp->next;
    }
    if (kount != clkp->nedges_total) return 6;
    return 0;
}

void clk_validate(CubicLatticeKnotPtr clkp, char* s)
{
    int code = clk_validate(clkp);
    if (code)
    {
        fprintf(stderr, "fails at %s with code %d\n", s, code);
        exit(102);
    }
}

bool clk_check_increment(ivector incr)
{
    // a valid increment should be +-1 in one direction only
    if (abs(incr[0]) + abs(incr[1]) + abs(incr[2]) == 1)
        return true;
    else
        return false;
}

bool clk_check_lattice(CubicLatticeKnotPtr knot)
{
    if (!knot) return false;
    // WARNING: this function has the side effect of clearing the knot path!
    //          use only at the end of a run to make sure things are OK

    ComponentCLKPtr comp = knot->fcomp;
    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        printf("checking lattice ... ");
        fflush(stdout);
        for (int e = 0; e < comp->nedges; e++)
        {
            if (knot->lattice[lat(ep->start)] != OCCUPIED)
            {
                printf("\n\nlattice location %d %d %d for edge with ID %d is not occupied!\n",
                       ep->start[0], ep->start[1], ep->start[2], ep->ID);
                fflush(stdout);
                exit(102);
            }
            knot->lattice[lat(ep->start)] = EMPTY;

            ep = ep->next;
        }

        if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
        {
            ivector endpoint;
            add_ivector(endpoint, comp->last_edge->start, comp->last_edge->increment);
            knot->lattice[lat(endpoint)] = EMPTY;
        }

        comp = comp->next;
    }

    int ix, iy, iz;
    ivector loc;

    // we've cleared the knot path, so every location in lattice should now be EMPTY

    for (ix = 0; ix < LATTICE_SIZE; ix++)
    {
        for (iy = 0; iy < LATTICE_SIZE; iy++)
        {
            for (iz = 0; iz < LATTICE_SIZE; iz++)
            {
                set_ivector(loc, ix, iy, iz);
                if (knot->lattice[lat(loc)] != EMPTY)
                {
                    printf("\n\nlattice location %d %d %d is not EMPTY as it should be!\n\n",
                           ix, iy, iz);
                    fflush(stdout);
                    exit(102);
                }
            }
        }
    }

    printf("lattice checks out!\n");
    return true;
}

void clk_check_increments(CubicLatticeKnotPtr knot)
{
    if (knot->nedges_total < 4)
    {
        fprintf(stderr,
                "\n\n*** clk_check_increments(): Ridiculous number of edges (%d) in knot!\n\n",
                knot->nedges_total);
        exit(102);
    }

    ComponentCLKPtr comp = knot->fcomp;
    FILE* fp = (FILE*)NULL;

    while (comp)
    {
        EdgePtr ep = comp->first_edge;
        ivector pos, start;

        copy_ivector(pos, comp->first_edge->start);
        copy_ivector(start, pos);

        int N = comp->nedges;
        if (comp->flags & COMPONENT_CLK_FLAG_OPEN) N--;

        for (int i = 0; i < N; i++)
        {
            if (!clk_check_increment(ep->increment))
            {
                fprintf(stderr, "\n\n clk_check_increments(): bad increment of (%d, %d, %d) seen!\n\n",
                        ep->increment[0], ep->increment[1], ep->increment[2]);
                exit(343);
            }
            if (clk_direction(ep->increment) != ep->dir)
            {
                fprintf(stderr, "\n\n clk_check_increments(): increment of (%d, %d, %d) is not direction %s!!\n\n",
                        ep->increment[0], ep->increment[1], ep->increment[2], clk_dir_name[ep->dir]);
                exit(93454);
            }
            add_ivector(pos, pos, ep->increment);
            if (!equal_ivector(pos, ep->next->start))
            {
                fprintf(stderr,
                        "\n\n *** clk_check_increments(): Bad increment!  Panicking!\n\n");
                exit(103);
            }
            if (ep->next)
            {
                if (ep->next->prev != ep)
                {
                    fprintf(stderr,
                            "\n\n *** clk_check_increments(): Edges not connected properly!  Panicking!\n\n");
                    exit(110);
                }
            }
            ep = ep->next;
        }
        if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN) && !equal_ivector(pos, start))
        {
            sub_ivector(pos, pos, knot->loffset);
            sub_ivector(start, start, knot->loffset);
            fprintf(stderr,
                    "\n\n *** clk_check_increments(): Component fails to close properly! (%d, %d, %d) != (%d, %d, %d) Panicking!\n\n",
                    start[0], start[1], start[2], pos[0], pos[1], pos[2]);

            clk_check_lattice(knot);
            exit(109);
        }
        comp = comp->next;
    }
}

bool clk_check_for_edge_hits(CubicLatticeKnotPtr knot, int dir, ivector test_location)
{
    // we should recentre the knot when an edge hit occurs
    // edge hits don't typically occur except for very small knots
    switch (dir)
    {
    case MOVE_NORTH:
        if (test_location[1] > LATTICE_SIZE - 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    case MOVE_SOUTH:
        if (test_location[1] < 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    case MOVE_EAST:
        if (test_location[0] > LATTICE_SIZE - 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    case MOVE_WEST:
        if (test_location[0] < 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    case MOVE_UP:
        if (test_location[2] > LATTICE_SIZE - 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    case MOVE_DOWN:
        if (test_location[2] < 2)
        {
            ++knot->edge_hits[dir];
            return true;
        }
        break;
    }

    return false;
}
