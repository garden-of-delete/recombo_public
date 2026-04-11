#include "writhe.h"

#include <cmath>

#ifdef WIN32
#define ISNAN _isnan
#else
#define ISNAN std::isnan
#endif

extern real rand_real(real low, real high);

static void zero_vector(vector3 v)
{
   v [0] = v [1] = v [2] = 0.0;
}

static void add_vector(vector3 v, vector3 w, vector3 z)
{
   v [0] = w [0] + z [0];
   v [1] = w [1] + z [1];
   v [2] = w [2] + z [2];
}

static void mult_vector(vector3 v, vector3 w, real con)
{
   v [0] = w [0] * con;
   v [1] = w [1] * con;
   v [2] = w [2] * con;
}

static real diff_distance_SQ(vector3 a, vector3 b)
{
   vector3 diff;
   sub_vector(diff, a, b);
   return diff [0] * diff [0] + diff [1] * diff [1] + diff [2] * diff [2];
}

double writhe(double &ACN, int nvert, ivector *vert, double jitter)
{
    vector3 V1, V2, V3, V4;
    vector3 R13, R14, R24, R23, R34, R12;
    vector3 N1, N2, N3, N4;
    vector3 R34xR12;

    double gauss_xing_number_sum = 0.0;
    double space_writhe_sum = 0.0;
    double exact_acn_contrib;

    int last = nvert - 1;

    for (int i = 0; i < nvert; i++)
    {
        set_vector(V1,
                (double) vert [i][0] + rand_real(-jitter, jitter),
                (double) vert [i][1] + rand_real(-jitter, jitter),
                (double) vert [i][2] + rand_real(-jitter, jitter));

        set_vector(V2,
                (double) vert [(i + 1) % nvert][0] + rand_real(-jitter, jitter),
                (double) vert [(i + 1) % nvert][1] + rand_real(-jitter, jitter),
                (double) vert [(i + 1) % nvert][2] + rand_real(-jitter, jitter));

        sub_vector(R12, V2, V1);

        for (int j = i + 2; j < last; j++)
        {
            set_vector(V3,
                    (double) vert [j][0] + rand_real(-jitter, jitter),
                    (double) vert [j][1] + rand_real(-jitter, jitter),
                    (double) vert [j][2] + rand_real(-jitter, jitter));

            set_vector(V4,
                    (double) vert [(j + 1) % nvert][0] + rand_real(-jitter, jitter),
                    (double) vert [(j + 1) % nvert][1] + rand_real(-jitter, jitter),
                    (double) vert [(j + 1) % nvert][2] + rand_real(-jitter, jitter));

            sub_vector(R13, V3, V1);
            sub_vector(R14, V4, V1);
            sub_vector(R24, V4, V2);
            sub_vector(R23, V3, V2);
            sub_vector(R34, V4, V3);

            cross(N1, R13, R14);
            normalize_vector(N1, N1);
            cross(N2, R14, R24);
            normalize_vector(N2, N2);
            cross(N3, R24, R23);
            normalize_vector(N3, N3);
            cross(N4, R23, R13);
            normalize_vector(N4, N4);

            exact_acn_contrib = asin(dot(N1, N2)) + asin(dot(N2, N3)) +
                                asin(dot(N3, N4)) + asin(dot(N4, N1));

            if (ISNAN(exact_acn_contrib))
                continue;

            gauss_xing_number_sum += exact_acn_contrib;
            cross(R34xR12, R34, R12);

            if (dot(R34xR12, R13) > 0.0)
                space_writhe_sum += exact_acn_contrib;
            else
                space_writhe_sum -= exact_acn_contrib;
        }
        last = nvert;
    }

    double space_writhe = space_writhe_sum / TWOPI;
    ACN = gauss_xing_number_sum / TWOPI;
    return space_writhe;
}

double writhe_open(ivector* AB, ivector* CD)
{
    double space_writhe = 0.0;
    vector3 A, B, C, D;
    vector3 R12, R13, R14, R23, R24, R34;
    vector3 N1, N2, N3, N4;
    vector3 R34xR12;

    set_vector(A, (double) AB[0][0], (double)AB[0][1], (double)AB[0][2]);
    set_vector(B, (double) AB[1][0], (double)AB[1][1], (double)AB[1][2]);
    set_vector(C, (double) CD[0][0], (double)CD[0][1], (double)CD[0][2]);
    set_vector(D, (double) CD[1][0], (double)CD[1][1], (double)CD[1][2]);

    sub_vector(R12, B, A);
    sub_vector(R13, C, A);
    sub_vector(R14, D, A);
    sub_vector(R23, C, B);
    sub_vector(R24, D, B);
    sub_vector(R34, D, C);

    cross(N1, R13, R14);
    cross(N2, R14, R24);
    cross(N3, R24, R23);
    cross(N4, R23, R13);

    normalize_vector(N1, N1);
    normalize_vector(N2, N2);
    normalize_vector(N3, N3);
    normalize_vector(N4, N4);

    double omega = asin(dot(N1, N2)) + asin(dot(N2, N3)) + asin(dot(N3, N4)) + asin(dot(N4, N1));
    cross(R34xR12, R34, R12);

    if (dot(R34xR12, R13) > 0.0)
        space_writhe += omega/TWOPI;
    else
        space_writhe -= omega/TWOPI;
    return space_writhe;
}

double radius_of_gyration(int nvert, ivector *vert)
{
   vector3 cofm, loc;
   zero_vector(cofm);
   for (int b = 0; b < nvert; b++)
   {
      set_vector(loc, (double) vert [b][0], (double) vert [b][1], (double) vert [b][2]);
      add_vector(cofm, cofm, loc);
   }
   mult_vector(cofm, cofm, 1.0 / (double) nvert);

   double sum = 0.0;
   for (int b = 0; b < nvert; b++)
   {
      set_vector(loc, (double) vert [b][0], (double) vert [b][1], (double) vert [b][2]);
      sum += diff_distance_SQ(loc, cofm);
   }

   return sqrt(sum / (double) nvert);
}
