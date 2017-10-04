#include "legacyBfacf.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

extern void set_sRand_seed_to_clocktime();
extern int rand_integer(int, int);
extern double rand_double(double low, double high);
extern double rand_uniform();
extern void sRandSimple(int seed);

// from findz.h
#define NUM_KNOTS 93

#define UNKNOWN_Z 108.0

#define NUM_KNOWN_LENGTHS 32
static int known_length [NUM_KNOWN_LENGTHS] = {56, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 124, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 300};

double zvalue_56 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 56
  0.196565332031,  // 3_1 good 56
  0.183241210937,  // 4_1 good 56
  0.173507421875,  // 5_1 good 56
  0.167922460937,  // 5_2 good 56
  0.147976171875,  // 6_1 good 56
  0.15021015625,  // 6_2 good 56
  0.146699609375,  // 6_3 good 56
  0.13106171875,  // 7_1 good 56
  0.12084921875,  // 7_2 good 56
  0.12212578125,  // 7_3 good 56
  0.1253171875,  // 7_4 good 56
  0.11191328125,  // 7_5 good 56
  0.1176578125,  // 7_6 good 56
  0.12340234375,  // 7_7 good 56
  0.07999921875,  // 8_1 good 56
  0.09021171875,  // 8_2 good 56
  0.083190625,  // 8_3 good 56
  0.09021171875,  // 8_4 good 56
  0.078084375,  // 8_5 good 56
  0.0844671875,  // 8_6 good 56
  0.09276484375,  // 8_7 good 56
  0.07999921875,  // 8_8 good 56
  0.098509375,  // 8_9 good 56
  0.078084375,  // 8_10 good 56
  0.078084375,  // 8_11 good 56
  0.0538296875,  // 8_12 good 56
  0.0946796875,  // 8_13 good 56
  0.0768078125,  // 8_14 good 56
  0.0742546875,  // 8_15 good 56
  0.083190625,  // 8_16 good 56
  0.078084375,  // 8_17 good 56
  0.083190625,  // 8_18 good 56
  0.14382734375,  // 8_19 good 56
  0.1278703125,  // 8_20 good 56
  0.1176578125,  // 8_21 good 56
  UNKNOWN_Z,  // 9_1 unknown 56
  UNKNOWN_Z,  // 9_2 unknown 56
  UNKNOWN_Z,  // 9_3 unknown 56
  UNKNOWN_Z,  // 9_4 unknown 56
  UNKNOWN_Z,  // 9_5 unknown 56
  UNKNOWN_Z,  // 9_6 unknown 56
  UNKNOWN_Z,  // 9_7 unknown 56
  UNKNOWN_Z,  // 9_8 unknown 56
  UNKNOWN_Z,  // 9_9 unknown 56
  UNKNOWN_Z,  // 9_10 unknown 56
  UNKNOWN_Z,  // 9_11 unknown 56
  UNKNOWN_Z,  // 9_12 unknown 56
  UNKNOWN_Z,  // 9_13 unknown 56
  UNKNOWN_Z,  // 9_14 unknown 56
  UNKNOWN_Z,  // 9_15 unknown 56
  UNKNOWN_Z,  // 9_16 unknown 56
  UNKNOWN_Z,  // 9_17 unknown 56
  UNKNOWN_Z,  // 9_18 unknown 56
  UNKNOWN_Z,  // 9_19 unknown 56
  UNKNOWN_Z,  // 9_20 unknown 56
  UNKNOWN_Z,  // 9_21 unknown 56
  UNKNOWN_Z,  // 9_22 unknown 56
  UNKNOWN_Z,  // 9_23 unknown 56
  UNKNOWN_Z,  // 9_24 unknown 56
  UNKNOWN_Z,  // 9_25 unknown 56
  UNKNOWN_Z,  // 9_26 unknown 56
  UNKNOWN_Z,  // 9_27 unknown 56
  UNKNOWN_Z,  // 9_28 unknown 56
  UNKNOWN_Z,  // 9_29 unknown 56
  UNKNOWN_Z,  // 9_30 unknown 56
  UNKNOWN_Z,  // 9_31 unknown 56
  UNKNOWN_Z,  // 9_32 unknown 56
  UNKNOWN_Z,  // 9_33 unknown 56
  UNKNOWN_Z,  // 9_34 unknown 56
  UNKNOWN_Z,  // 9_35 unknown 56
  UNKNOWN_Z,  // 9_36 unknown 56
  UNKNOWN_Z,  // 9_37 unknown 56
  UNKNOWN_Z,  // 9_38 unknown 56
  UNKNOWN_Z,  // 9_39 unknown 56
  UNKNOWN_Z,  // 9_40 unknown 56
  UNKNOWN_Z,  // 9_41 unknown 56
  UNKNOWN_Z,  // 9_42 unknown 56
  UNKNOWN_Z,  // 9_43 unknown 56
  UNKNOWN_Z,  // 9_44 unknown 56
  UNKNOWN_Z,  // 9_45 unknown 56
  UNKNOWN_Z,  // 9_46 unknown 56
  UNKNOWN_Z,  // 9_47 unknown 56
  UNKNOWN_Z,  // 9_48 unknown 56
  UNKNOWN_Z,  // 9_49 unknown 56
  0.15021015625,  // 6c_1 good 56
  0.152444140625,  // 6c_2 good 56
  0.1176578125,  // 7c_1 good 56
  0.0870203125,  // 8c_1 good 56
  0.0921265625,  // 8c_2 good 56
  0.0742546875,  // 8c_3 good 56
  0.078084375,  // 8c_4 good 56
  0.0691484375}; // 8c_5 good 56

double zvalue_60 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 60
  0.198787988281,  // 3_1 good 60
  0.187624902344,  // 4_1 good 60
  0.179172851562,  // 5_1 good 60
  0.174867089844,  // 5_2 good 60
  0.159079296875,  // 6_1 good 60
  0.159876660156,  // 6_2 good 60
  0.157325097656,  // 6_3 good 60
  0.14568359375,  // 7_1 good 60
  0.13802890625,  // 7_2 good 60
  0.13802890625,  // 7_3 good 60
  0.139942578125,  // 7_4 good 60
  0.131331054688,  // 7_5 good 60
  0.133563671875,  // 7_6 good 60
  0.13802890625,  // 7_7 good 60
  0.10996171875,  // 8_1 good 60
  0.11506484375,  // 8_2 good 60
  0.109323828125,  // 8_3 good 60
  0.11506484375,  // 8_4 good 60
  0.10741015625,  // 8_5 good 60
  0.10741015625,  // 8_6 good 60
  0.11506484375,  // 8_7 good 60
  0.10485859375,  // 8_8 good 60
  0.118254296875,  // 8_9 good 60
  0.1035828125,  // 8_10 good 60
  0.10485859375,  // 8_11 good 60
  0.090187109375,  // 8_12 good 60
  0.11506484375,  // 8_13 good 60
  0.102944921875,  // 8_14 good 60
  0.101669140625,  // 8_15 good 60
  0.111875390625,  // 8_16 good 60
  0.10741015625,  // 8_17 good 60
  0.110599609375,  // 8_18 good 60
  0.15588984375,  // 8_19 good 60
  0.143769921875,  // 8_20 good 60
  0.13547734375,  // 8_21 good 60
  UNKNOWN_Z,  // 9_1 unknown 60
  UNKNOWN_Z,  // 9_2 unknown 60
  UNKNOWN_Z,  // 9_3 unknown 60
  UNKNOWN_Z,  // 9_4 unknown 60
  UNKNOWN_Z,  // 9_5 unknown 60
  UNKNOWN_Z,  // 9_6 unknown 60
  UNKNOWN_Z,  // 9_7 unknown 60
  UNKNOWN_Z,  // 9_8 unknown 60
  UNKNOWN_Z,  // 9_9 unknown 60
  UNKNOWN_Z,  // 9_10 unknown 60
  UNKNOWN_Z,  // 9_11 unknown 60
  UNKNOWN_Z,  // 9_12 unknown 60
  UNKNOWN_Z,  // 9_13 unknown 60
  UNKNOWN_Z,  // 9_14 unknown 60
  UNKNOWN_Z,  // 9_15 unknown 60
  UNKNOWN_Z,  // 9_16 unknown 60
  UNKNOWN_Z,  // 9_17 unknown 60
  UNKNOWN_Z,  // 9_18 unknown 60
  UNKNOWN_Z,  // 9_19 unknown 60
  UNKNOWN_Z,  // 9_20 unknown 60
  UNKNOWN_Z,  // 9_21 unknown 60
  UNKNOWN_Z,  // 9_22 unknown 60
  UNKNOWN_Z,  // 9_23 unknown 60
  UNKNOWN_Z,  // 9_24 unknown 60
  UNKNOWN_Z,  // 9_25 unknown 60
  UNKNOWN_Z,  // 9_26 unknown 60
  UNKNOWN_Z,  // 9_27 unknown 60
  UNKNOWN_Z,  // 9_28 unknown 60
  UNKNOWN_Z,  // 9_29 unknown 60
  UNKNOWN_Z,  // 9_30 unknown 60
  UNKNOWN_Z,  // 9_31 unknown 60
  UNKNOWN_Z,  // 9_32 unknown 60
  UNKNOWN_Z,  // 9_33 unknown 60
  UNKNOWN_Z,  // 9_34 unknown 60
  UNKNOWN_Z,  // 9_35 unknown 60
  UNKNOWN_Z,  // 9_36 unknown 60
  UNKNOWN_Z,  // 9_37 unknown 60
  UNKNOWN_Z,  // 9_38 unknown 60
  UNKNOWN_Z,  // 9_39 unknown 60
  UNKNOWN_Z,  // 9_40 unknown 60
  UNKNOWN_Z,  // 9_41 unknown 60
  UNKNOWN_Z,  // 9_42 unknown 60
  UNKNOWN_Z,  // 9_43 unknown 60
  UNKNOWN_Z,  // 9_44 unknown 60
  UNKNOWN_Z,  // 9_45 unknown 60
  UNKNOWN_Z,  // 9_46 unknown 60
  UNKNOWN_Z,  // 9_47 unknown 60
  UNKNOWN_Z,  // 9_48 unknown 60
  UNKNOWN_Z,  // 9_49 unknown 60
  0.160355078125,  // 6c_1 good 60
  0.161949804687,  // 6c_2 good 60
  0.13547734375,  // 7c_1 good 60
  0.113151171875,  // 8c_1 good 60
  0.11761640625,  // 8c_2 good 60
  0.10485859375,  // 8c_3 good 60
  0.10485859375,  // 8c_4 good 60
  0.0984796875}; // 8c_5 good 60

double zvalue_65 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 65
  0.201020605469,  // 3_1 good 65
  0.191691455078,  // 4_1 good 65
  0.185073339844,  // 5_1 good 65
  0.181724414062,  // 5_2 good 65
  0.168966601562,  // 6_1 good 65
  0.168966601562,  // 6_2 good 65
  0.167052929688,  // 6_3 good 65
  0.159079296875,  // 7_1 good 65
  0.152381445312,  // 7_2 good 65
  0.153019335938,  // 7_3 good 65
  0.153976171875,  // 7_4 good 65
  0.147916210938,  // 7_5 good 65
  0.14823515625,  // 7_6 good 65
  0.151105664062,  // 7_7 good 65
  0.131968945312,  // 8_1 good 65
  0.134839453125,  // 8_2 good 65
  0.131968945312,  // 8_3 good 65
  0.134839453125,  // 8_4 good 65
  0.129736328125,  // 8_5 good 65
  0.12782265625,  // 8_6 good 65
  0.133563671875,  // 8_7 good 65
  0.12782265625,  // 8_8 good 65
  0.13547734375,  // 8_9 good 65
  0.12527109375,  // 8_10 good 65
  0.12782265625,  // 8_11 good 65
  0.11761640625,  // 8_12 good 65
  0.133563671875,  // 8_13 good 65
  0.12527109375,  // 8_14 good 65
  0.12527109375,  // 8_15 good 65
  0.133563671875,  // 8_16 good 65
  0.129736328125,  // 8_17 good 65
  0.131968945312,  // 8_18 good 65
  0.166415039062,  // 8_19 good 65
  0.158122460937,  // 8_20 good 65
  0.151105664062,  // 8_21 good 65
  UNKNOWN_Z,  // 9_1 unknown 65
  UNKNOWN_Z,  // 9_2 unknown 65
  UNKNOWN_Z,  // 9_3 unknown 65
  UNKNOWN_Z,  // 9_4 unknown 65
  UNKNOWN_Z,  // 9_5 unknown 65
  UNKNOWN_Z,  // 9_6 unknown 65
  UNKNOWN_Z,  // 9_7 unknown 65
  UNKNOWN_Z,  // 9_8 unknown 65
  UNKNOWN_Z,  // 9_9 unknown 65
  UNKNOWN_Z,  // 9_10 unknown 65
  UNKNOWN_Z,  // 9_11 unknown 65
  UNKNOWN_Z,  // 9_12 unknown 65
  UNKNOWN_Z,  // 9_13 unknown 65
  UNKNOWN_Z,  // 9_14 unknown 65
  UNKNOWN_Z,  // 9_15 unknown 65
  UNKNOWN_Z,  // 9_16 unknown 65
  UNKNOWN_Z,  // 9_17 unknown 65
  UNKNOWN_Z,  // 9_18 unknown 65
  UNKNOWN_Z,  // 9_19 unknown 65
  UNKNOWN_Z,  // 9_20 unknown 65
  UNKNOWN_Z,  // 9_21 unknown 65
  UNKNOWN_Z,  // 9_22 unknown 65
  UNKNOWN_Z,  // 9_23 unknown 65
  UNKNOWN_Z,  // 9_24 unknown 65
  UNKNOWN_Z,  // 9_25 unknown 65
  UNKNOWN_Z,  // 9_26 unknown 65
  UNKNOWN_Z,  // 9_27 unknown 65
  UNKNOWN_Z,  // 9_28 unknown 65
  UNKNOWN_Z,  // 9_29 unknown 65
  UNKNOWN_Z,  // 9_30 unknown 65
  UNKNOWN_Z,  // 9_31 unknown 65
  UNKNOWN_Z,  // 9_32 unknown 65
  UNKNOWN_Z,  // 9_33 unknown 65
  UNKNOWN_Z,  // 9_34 unknown 65
  UNKNOWN_Z,  // 9_35 unknown 65
  UNKNOWN_Z,  // 9_36 unknown 65
  UNKNOWN_Z,  // 9_37 unknown 65
  UNKNOWN_Z,  // 9_38 unknown 65
  UNKNOWN_Z,  // 9_39 unknown 65
  UNKNOWN_Z,  // 9_40 unknown 65
  UNKNOWN_Z,  // 9_41 unknown 65
  UNKNOWN_Z,  // 9_42 unknown 65
  UNKNOWN_Z,  // 9_43 unknown 65
  UNKNOWN_Z,  // 9_44 unknown 65
  UNKNOWN_Z,  // 9_45 unknown 65
  UNKNOWN_Z,  // 9_46 unknown 65
  UNKNOWN_Z,  // 9_47 unknown 65
  UNKNOWN_Z,  // 9_48 unknown 65
  UNKNOWN_Z,  // 9_49 unknown 65
  0.170082910156,  // 6c_1 good 65
  0.170561328125,  // 6c_2 good 65
  0.150148828125,  // 7c_1 good 65
  0.134839453125,  // 8c_1 good 65
  0.136434179687,  // 8c_2 good 65
  0.12782265625,  // 8c_3 good 65
  0.12782265625,  // 8c_4 good 65
  0.122081640625}; // 8c_5 good 65

double zvalue_70 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 70
  0.202695068359,  // 3_1 good 70
  0.194960644531,  // 4_1 good 70
  0.189379101562,  // 5_1 good 70
  0.186189648438,  // 5_2 good 70
  0.175983398437,  // 6_1 good 70
  0.175983398437,  // 6_2 good 70
  0.174388671875,  // 6_3 good 70
  0.168966601562,  // 7_1 good 70
  0.163225585938,  // 7_2 good 70
  0.163225585938,  // 7_3 good 70
  0.163225585938,  // 7_4 good 70
  0.159079296875,  // 7_5 good 70
  0.159079296875,  // 7_6 good 70
  0.160674023437,  // 7_7 good 70
  0.147278320312,  // 8_1 good 70
  0.14823515625,  // 8_2 good 70
  0.146640429688,  // 8_3 good 70
  0.14823515625,  // 8_4 good 70
  0.145045703125,  // 8_5 good 70
  0.142494140625,  // 8_6 good 70
  0.146640429688,  // 8_7 good 70
  0.142494140625,  // 8_8 good 70
  0.147916210938,  // 8_9 good 70
  0.141537304687,  // 8_10 good 70
  0.143769921875,  // 8_11 good 70
  0.137072070312,  // 8_12 good 70
  0.146640429688,  // 8_13 good 70
  0.141218359375,  // 8_14 good 70
  0.139942578125,  // 8_15 good 70
  0.147916210938,  // 8_16 good 70
  0.145045703125,  // 8_17 good 70
  0.146002539062,  // 8_18 good 70
  0.174548144531,  // 8_19 good 70
  0.167531347656,  // 8_20 good 70
  0.161949804687,  // 8_21 good 70
  UNKNOWN_Z,  // 9_1 unknown 70
  UNKNOWN_Z,  // 9_2 unknown 70
  UNKNOWN_Z,  // 9_3 unknown 70
  UNKNOWN_Z,  // 9_4 unknown 70
  UNKNOWN_Z,  // 9_5 unknown 70
  UNKNOWN_Z,  // 9_6 unknown 70
  UNKNOWN_Z,  // 9_7 unknown 70
  UNKNOWN_Z,  // 9_8 unknown 70
  UNKNOWN_Z,  // 9_9 unknown 70
  UNKNOWN_Z,  // 9_10 unknown 70
  UNKNOWN_Z,  // 9_11 unknown 70
  UNKNOWN_Z,  // 9_12 unknown 70
  UNKNOWN_Z,  // 9_13 unknown 70
  UNKNOWN_Z,  // 9_14 unknown 70
  UNKNOWN_Z,  // 9_15 unknown 70
  UNKNOWN_Z,  // 9_16 unknown 70
  UNKNOWN_Z,  // 9_17 unknown 70
  UNKNOWN_Z,  // 9_18 unknown 70
  UNKNOWN_Z,  // 9_19 unknown 70
  UNKNOWN_Z,  // 9_20 unknown 70
  UNKNOWN_Z,  // 9_21 unknown 70
  UNKNOWN_Z,  // 9_22 unknown 70
  UNKNOWN_Z,  // 9_23 unknown 70
  UNKNOWN_Z,  // 9_24 unknown 70
  UNKNOWN_Z,  // 9_25 unknown 70
  UNKNOWN_Z,  // 9_26 unknown 70
  UNKNOWN_Z,  // 9_27 unknown 70
  UNKNOWN_Z,  // 9_28 unknown 70
  UNKNOWN_Z,  // 9_29 unknown 70
  UNKNOWN_Z,  // 9_30 unknown 70
  UNKNOWN_Z,  // 9_31 unknown 70
  UNKNOWN_Z,  // 9_32 unknown 70
  UNKNOWN_Z,  // 9_33 unknown 70
  UNKNOWN_Z,  // 9_34 unknown 70
  UNKNOWN_Z,  // 9_35 unknown 70
  UNKNOWN_Z,  // 9_36 unknown 70
  UNKNOWN_Z,  // 9_37 unknown 70
  UNKNOWN_Z,  // 9_38 unknown 70
  UNKNOWN_Z,  // 9_39 unknown 70
  UNKNOWN_Z,  // 9_40 unknown 70
  UNKNOWN_Z,  // 9_41 unknown 70
  UNKNOWN_Z,  // 9_42 unknown 70
  UNKNOWN_Z,  // 9_43 unknown 70
  UNKNOWN_Z,  // 9_44 unknown 70
  UNKNOWN_Z,  // 9_45 unknown 70
  UNKNOWN_Z,  // 9_46 unknown 70
  UNKNOWN_Z,  // 9_47 unknown 70
  UNKNOWN_Z,  // 9_48 unknown 70
  UNKNOWN_Z,  // 9_49 unknown 70
  0.176621289062,  // 6c_1 good 70
  0.177418652344,  // 6c_2 good 70
  0.161311914062,  // 7c_1 good 70
  0.148873046875,  // 8c_1 good 70
  0.150148828125,  // 8c_2 good 70
  0.143769921875,  // 8c_3 good 70
  0.143769921875,  // 8c_4 good 70
  0.139942578125}; // 8c_5 good 70

double zvalue_75 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 75
  0.204064941406,  // 3_1 good 75
  0.197488708496,  // 4_1 good 75
  0.192731933594,  // 5_1 good 75
  0.190261230469,  // 5_2 good 75
  0.181721191406,  // 6_1 good 75
  0.181153869629,  // 6_2 good 75
  0.179989013672,  // 6_3 good 75
  0.175598144531,  // 7_1 good 75
  0.171274414062,  // 7_2 good 75
  0.170871582031,  // 7_3 good 75
  0.170737304687,  // 7_4 good 75
  0.167836914063,  // 7_5 good 75
  0.167622070313,  // 7_6 good 75
  0.168548583984,  // 7_7 good 75
  0.158383789063,  // 8_1 good 75
  0.158598632813,  // 8_2 good 75
  0.157846679687,  // 8_3 good 75
  0.158491210938,  // 8_4 good 75
  0.155913085938,  // 8_5 good 75
  0.154194335937,  // 8_6 good 75
  0.156772460938,  // 8_7 good 75
  0.154194335937,  // 8_8 good 75
  0.157416992187,  // 8_9 good 75
  0.153442382812,  // 8_10 good 75
  0.154409179687,  // 8_11 good 75
  0.149897460938,  // 8_12 good 75
  0.156396484375,  // 8_13 good 75
  0.152475585937,  // 8_14 good 75
  0.152099609375,  // 8_15 good 75
  0.157900390625,  // 8_16 good 75
  0.155698242188,  // 8_17 good 75
  0.156665039063,  // 8_18 good 75
  0.180512695313,  // 8_19 good 75
  0.174685058594,  // 8_20 good 75
  0.169931640625,  // 8_21 good 75
  UNKNOWN_Z,  // 9_1 unknown 75
  UNKNOWN_Z,  // 9_2 unknown 75
  UNKNOWN_Z,  // 9_3 unknown 75
  UNKNOWN_Z,  // 9_4 unknown 75
  UNKNOWN_Z,  // 9_5 unknown 75
  UNKNOWN_Z,  // 9_6 unknown 75
  UNKNOWN_Z,  // 9_7 unknown 75
  UNKNOWN_Z,  // 9_8 unknown 75
  UNKNOWN_Z,  // 9_9 unknown 75
  UNKNOWN_Z,  // 9_10 unknown 75
  UNKNOWN_Z,  // 9_11 unknown 75
  UNKNOWN_Z,  // 9_12 unknown 75
  UNKNOWN_Z,  // 9_13 unknown 75
  UNKNOWN_Z,  // 9_14 unknown 75
  UNKNOWN_Z,  // 9_15 unknown 75
  UNKNOWN_Z,  // 9_16 unknown 75
  UNKNOWN_Z,  // 9_17 unknown 75
  UNKNOWN_Z,  // 9_18 unknown 75
  UNKNOWN_Z,  // 9_19 unknown 75
  UNKNOWN_Z,  // 9_20 unknown 75
  UNKNOWN_Z,  // 9_21 unknown 75
  UNKNOWN_Z,  // 9_22 unknown 75
  UNKNOWN_Z,  // 9_23 unknown 75
  UNKNOWN_Z,  // 9_24 unknown 75
  UNKNOWN_Z,  // 9_25 unknown 75
  UNKNOWN_Z,  // 9_26 unknown 75
  UNKNOWN_Z,  // 9_27 unknown 75
  UNKNOWN_Z,  // 9_28 unknown 75
  UNKNOWN_Z,  // 9_29 unknown 75
  UNKNOWN_Z,  // 9_30 unknown 75
  UNKNOWN_Z,  // 9_31 unknown 75
  UNKNOWN_Z,  // 9_32 unknown 75
  UNKNOWN_Z,  // 9_33 unknown 75
  UNKNOWN_Z,  // 9_34 unknown 75
  UNKNOWN_Z,  // 9_35 unknown 75
  UNKNOWN_Z,  // 9_36 unknown 75
  UNKNOWN_Z,  // 9_37 unknown 75
  UNKNOWN_Z,  // 9_38 unknown 75
  UNKNOWN_Z,  // 9_39 unknown 75
  UNKNOWN_Z,  // 9_40 unknown 75
  UNKNOWN_Z,  // 9_41 unknown 75
  UNKNOWN_Z,  // 9_42 unknown 75
  UNKNOWN_Z,  // 9_43 unknown 75
  UNKNOWN_Z,  // 9_44 unknown 75
  UNKNOWN_Z,  // 9_45 unknown 75
  UNKNOWN_Z,  // 9_46 unknown 75
  UNKNOWN_Z,  // 9_47 unknown 75
  UNKNOWN_Z,  // 9_48 unknown 75
  UNKNOWN_Z,  // 9_49 unknown 75
  0.182283789063,  // 6c_1 good 75
  0.182283789063,  // 6c_2 good 75
  0.169039453125,  // 7c_1 good 75
  0.15914609375,  // 8c_1 good 75
  0.159943945313,  // 8c_2 good 75
  0.15531640625,  // 8c_3 good 75
  0.15531640625,  // 8c_4 good 75
  0.151805859375}; // 8c_5 good 75

double zvalue_80 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 80
  0.205206762695,  // 3_1 good 80
  0.199505615234,  // 4_1 good 80
  0.195359326172,  // 5_1 good 80
  0.193365917969,  // 5_2 good 80
  0.186189648438,  // 6_1 good 80
  0.185392285156,  // 6_2 good 80
  0.184275976562,  // 6_3 good 80
  0.181086523437,  // 7_1 good 80
  0.177099707031,  // 7_2 good 80
  0.176621289062,  // 7_3 good 80
  0.176621289062,  // 7_4 good 80
  0.174388671875,  // 7_5 good 80
  0.174069726562,  // 7_6 good 80
  0.174548144531,  // 7_7 good 80
  0.166415039062,  // 8_1 good 80
  0.166415039062,  // 8_2 good 80
  0.166415039062,  // 8_3 good 80
  0.166415039062,  // 8_4 good 80
  0.164182421875,  // 8_5 good 80
  0.163225585938,  // 8_6 good 80
  0.164660839844,  // 8_7 good 80
  0.162587695313,  // 8_8 good 80
  0.164979785156,  // 8_9 good 80
  0.162428222656,  // 8_10 good 80
  0.162587695313,  // 8_11 good 80
  0.160355078125,  // 8_12 good 80
  0.164182421875,  // 8_13 good 80
  0.161311914062,  // 8_14 good 80
  0.161311914062,  // 8_15 good 80
  0.165458203125,  // 8_16 good 80
  0.163863476562,  // 8_17 good 80
  0.164660839844,  // 8_18 good 80
  0.185073339844,  // 8_19 good 80
  0.180289160156,  // 8_20 good 80
  0.175983398437,  // 8_21 good 80
  UNKNOWN_Z,  // 9_1 unknown 80
  UNKNOWN_Z,  // 9_2 unknown 80
  UNKNOWN_Z,  // 9_3 unknown 80
  UNKNOWN_Z,  // 9_4 unknown 80
  UNKNOWN_Z,  // 9_5 unknown 80
  UNKNOWN_Z,  // 9_6 unknown 80
  UNKNOWN_Z,  // 9_7 unknown 80
  UNKNOWN_Z,  // 9_8 unknown 80
  UNKNOWN_Z,  // 9_9 unknown 80
  UNKNOWN_Z,  // 9_10 unknown 80
  UNKNOWN_Z,  // 9_11 unknown 80
  UNKNOWN_Z,  // 9_12 unknown 80
  UNKNOWN_Z,  // 9_13 unknown 80
  UNKNOWN_Z,  // 9_14 unknown 80
  UNKNOWN_Z,  // 9_15 unknown 80
  UNKNOWN_Z,  // 9_16 unknown 80
  UNKNOWN_Z,  // 9_17 unknown 80
  UNKNOWN_Z,  // 9_18 unknown 80
  UNKNOWN_Z,  // 9_19 unknown 80
  UNKNOWN_Z,  // 9_20 unknown 80
  UNKNOWN_Z,  // 9_21 unknown 80
  UNKNOWN_Z,  // 9_22 unknown 80
  UNKNOWN_Z,  // 9_23 unknown 80
  UNKNOWN_Z,  // 9_24 unknown 80
  UNKNOWN_Z,  // 9_25 unknown 80
  UNKNOWN_Z,  // 9_26 unknown 80
  UNKNOWN_Z,  // 9_27 unknown 80
  UNKNOWN_Z,  // 9_28 unknown 80
  UNKNOWN_Z,  // 9_29 unknown 80
  UNKNOWN_Z,  // 9_30 unknown 80
  UNKNOWN_Z,  // 9_31 unknown 80
  UNKNOWN_Z,  // 9_32 unknown 80
  UNKNOWN_Z,  // 9_33 unknown 80
  UNKNOWN_Z,  // 9_34 unknown 80
  UNKNOWN_Z,  // 9_35 unknown 80
  UNKNOWN_Z,  // 9_36 unknown 80
  UNKNOWN_Z,  // 9_37 unknown 80
  UNKNOWN_Z,  // 9_38 unknown 80
  UNKNOWN_Z,  // 9_39 unknown 80
  UNKNOWN_Z,  // 9_40 unknown 80
  UNKNOWN_Z,  // 9_41 unknown 80
  UNKNOWN_Z,  // 9_42 unknown 80
  UNKNOWN_Z,  // 9_43 unknown 80
  UNKNOWN_Z,  // 9_44 unknown 80
  UNKNOWN_Z,  // 9_45 unknown 80
  UNKNOWN_Z,  // 9_46 unknown 80
  UNKNOWN_Z,  // 9_47 unknown 80
  UNKNOWN_Z,  // 9_48 unknown 80
  UNKNOWN_Z,  // 9_49 unknown 80
  0.186189648438,  // 6c_1 good 80
  0.186189648438,  // 6c_2 good 80
  0.175186035156,  // 7c_1 good 80
  0.167052929688,  // 8c_1 good 80
  0.167531347656,  // 8c_2 good 80
  0.164182421875,  // 8c_3 good 80
  0.164182421875,  // 8c_4 good 80
  0.161311914062}; // 8c_5 good 80

double zvalue_85 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 85
  0.206083862305,  // 3_1 good 85
  0.201140209961,  // 4_1 good 85
  0.197512207031,  // 5_1 good 85
  0.195678271484,  // 5_2 good 85
  0.189379101562,  // 6_1 good 85
  0.188820947266,  // 6_2 good 85
  0.187943847656,  // 6_3 good 85
  0.185073339844,  // 7_1 good 85
  0.181724414062,  // 7_2 good 85
  0.181724414062,  // 7_3 good 85
  0.181086523437,  // 7_4 good 85
  0.179491796875,  // 7_5 good 85
  0.179172851562,  // 7_6 good 85
  0.179172851562,  // 7_7 good 85
  0.172793945312,  // 8_1 good 85
  0.172634472656,  // 8_2 good 85
  0.172634472656,  // 8_3 good 85
  0.172315527344,  // 8_4 good 85
  0.170880273438,  // 8_5 good 85
  0.169763964844,  // 8_6 good 85
  0.170880273438,  // 8_7 good 85
  0.169445019531,  // 8_8 good 85
  0.170880273438,  // 8_9 good 85
  0.169285546875,  // 8_10 good 85
  0.169285546875,  // 8_11 good 85
  0.167212402344,  // 8_12 good 85
  0.170401855469,  // 8_13 good 85
  0.168328710938,  // 8_14 good 85
  0.167690820312,  // 8_15 good 85
  0.171518164062,  // 8_16 good 85
  0.170082910156,  // 8_17 good 85
  0.170561328125,  // 8_18 good 85
  0.188741210938,  // 8_19 good 85
  0.184275976562,  // 8_20 good 85
  0.181086523437,  // 8_21 good 85
  UNKNOWN_Z,  // 9_1 unknown 85
  UNKNOWN_Z,  // 9_2 unknown 85
  UNKNOWN_Z,  // 9_3 unknown 85
  UNKNOWN_Z,  // 9_4 unknown 85
  UNKNOWN_Z,  // 9_5 unknown 85
  UNKNOWN_Z,  // 9_6 unknown 85
  UNKNOWN_Z,  // 9_7 unknown 85
  UNKNOWN_Z,  // 9_8 unknown 85
  UNKNOWN_Z,  // 9_9 unknown 85
  UNKNOWN_Z,  // 9_10 unknown 85
  UNKNOWN_Z,  // 9_11 unknown 85
  UNKNOWN_Z,  // 9_12 unknown 85
  UNKNOWN_Z,  // 9_13 unknown 85
  UNKNOWN_Z,  // 9_14 unknown 85
  UNKNOWN_Z,  // 9_15 unknown 85
  UNKNOWN_Z,  // 9_16 unknown 85
  UNKNOWN_Z,  // 9_17 unknown 85
  UNKNOWN_Z,  // 9_18 unknown 85
  UNKNOWN_Z,  // 9_19 unknown 85
  UNKNOWN_Z,  // 9_20 unknown 85
  UNKNOWN_Z,  // 9_21 unknown 85
  UNKNOWN_Z,  // 9_22 unknown 85
  UNKNOWN_Z,  // 9_23 unknown 85
  UNKNOWN_Z,  // 9_24 unknown 85
  UNKNOWN_Z,  // 9_25 unknown 85
  UNKNOWN_Z,  // 9_26 unknown 85
  UNKNOWN_Z,  // 9_27 unknown 85
  UNKNOWN_Z,  // 9_28 unknown 85
  UNKNOWN_Z,  // 9_29 unknown 85
  UNKNOWN_Z,  // 9_30 unknown 85
  UNKNOWN_Z,  // 9_31 unknown 85
  UNKNOWN_Z,  // 9_32 unknown 85
  UNKNOWN_Z,  // 9_33 unknown 85
  UNKNOWN_Z,  // 9_34 unknown 85
  UNKNOWN_Z,  // 9_35 unknown 85
  UNKNOWN_Z,  // 9_36 unknown 85
  UNKNOWN_Z,  // 9_37 unknown 85
  UNKNOWN_Z,  // 9_38 unknown 85
  UNKNOWN_Z,  // 9_39 unknown 85
  UNKNOWN_Z,  // 9_40 unknown 85
  UNKNOWN_Z,  // 9_41 unknown 85
  UNKNOWN_Z,  // 9_42 unknown 85
  UNKNOWN_Z,  // 9_43 unknown 85
  UNKNOWN_Z,  // 9_44 unknown 85
  UNKNOWN_Z,  // 9_45 unknown 85
  UNKNOWN_Z,  // 9_46 unknown 85
  UNKNOWN_Z,  // 9_47 unknown 85
  UNKNOWN_Z,  // 9_48 unknown 85
  UNKNOWN_Z,  // 9_49 unknown 85
  0.189379101562,  // 6c_1 good 85
  0.189538574219,  // 6c_2 good 85
  0.179970214844,  // 7c_1 good 85
  0.173431835938,  // 8c_1 good 85
  0.173431835938,  // 8c_2 good 85
  0.170561328125,  // 8c_3 good 85
  0.170561328125,  // 8c_4 good 85
  0.168328710938}; // 8c_5 good 85

double zvalue_90 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 90
  0.206881225586,  // 3_1 good 90
  0.202415991211,  // 4_1 good 90
  0.199306274414,  // 5_1 good 90
  0.197751416016,  // 5_2 good 90
  0.192409082031,  // 6_1 good 90
  0.191691455078,  // 6_2 good 90
  0.190814355469,  // 6_3 good 90
  0.188581738281,  // 7_1 good 90
  0.185711230469,  // 7_2 good 90
  0.185392285156,  // 7_3 good 90
  0.185073339844,  // 7_4 good 90
  0.183638085938,  // 7_5 good 90
  0.183159667969,  // 7_6 good 90
  0.183478613281,  // 7_7 good 90
  0.177897070312,  // 8_1 good 90
  0.177737597656,  // 8_2 good 90
  0.177418652344,  // 8_3 good 90
  0.177418652344,  // 8_4 good 90
  0.175983398437,  // 8_5 good 90
  0.175186035156,  // 8_6 good 90
  0.176621289062,  // 8_7 good 90
  0.174867089844,  // 8_8 good 90
  0.175983398437,  // 8_9 good 90
  0.174548144531,  // 8_10 good 90
  0.174867089844,  // 8_11 good 90
  0.173431835938,  // 8_12 good 90
  0.175504980469,  // 8_13 good 90
  0.174069726562,  // 8_14 good 90
  0.173431835938,  // 8_15 good 90
  0.175983398437,  // 8_16 good 90
  0.175186035156,  // 8_17 good 90
  0.175664453125,  // 8_18 good 90
  0.191531982422,  // 8_19 good 90
  0.187943847656,  // 8_20 good 90
  0.185073339844,  // 8_21 good 90
  UNKNOWN_Z,  // 9_1 unknown 90
  UNKNOWN_Z,  // 9_2 unknown 90
  UNKNOWN_Z,  // 9_3 unknown 90
  UNKNOWN_Z,  // 9_4 unknown 90
  UNKNOWN_Z,  // 9_5 unknown 90
  UNKNOWN_Z,  // 9_6 unknown 90
  UNKNOWN_Z,  // 9_7 unknown 90
  UNKNOWN_Z,  // 9_8 unknown 90
  UNKNOWN_Z,  // 9_9 unknown 90
  UNKNOWN_Z,  // 9_10 unknown 90
  UNKNOWN_Z,  // 9_11 unknown 90
  UNKNOWN_Z,  // 9_12 unknown 90
  UNKNOWN_Z,  // 9_13 unknown 90
  UNKNOWN_Z,  // 9_14 unknown 90
  UNKNOWN_Z,  // 9_15 unknown 90
  UNKNOWN_Z,  // 9_16 unknown 90
  UNKNOWN_Z,  // 9_17 unknown 90
  UNKNOWN_Z,  // 9_18 unknown 90
  UNKNOWN_Z,  // 9_19 unknown 90
  UNKNOWN_Z,  // 9_20 unknown 90
  UNKNOWN_Z,  // 9_21 unknown 90
  UNKNOWN_Z,  // 9_22 unknown 90
  UNKNOWN_Z,  // 9_23 unknown 90
  UNKNOWN_Z,  // 9_24 unknown 90
  UNKNOWN_Z,  // 9_25 unknown 90
  UNKNOWN_Z,  // 9_26 unknown 90
  UNKNOWN_Z,  // 9_27 unknown 90
  UNKNOWN_Z,  // 9_28 unknown 90
  UNKNOWN_Z,  // 9_29 unknown 90
  UNKNOWN_Z,  // 9_30 unknown 90
  UNKNOWN_Z,  // 9_31 unknown 90
  UNKNOWN_Z,  // 9_32 unknown 90
  UNKNOWN_Z,  // 9_33 unknown 90
  UNKNOWN_Z,  // 9_34 unknown 90
  UNKNOWN_Z,  // 9_35 unknown 90
  UNKNOWN_Z,  // 9_36 unknown 90
  UNKNOWN_Z,  // 9_37 unknown 90
  UNKNOWN_Z,  // 9_38 unknown 90
  UNKNOWN_Z,  // 9_39 unknown 90
  UNKNOWN_Z,  // 9_40 unknown 90
  UNKNOWN_Z,  // 9_41 unknown 90
  UNKNOWN_Z,  // 9_42 unknown 90
  UNKNOWN_Z,  // 9_43 unknown 90
  UNKNOWN_Z,  // 9_44 unknown 90
  UNKNOWN_Z,  // 9_45 unknown 90
  UNKNOWN_Z,  // 9_46 unknown 90
  UNKNOWN_Z,  // 9_47 unknown 90
  UNKNOWN_Z,  // 9_48 unknown 90
  UNKNOWN_Z,  // 9_49 unknown 90
  0.192090136719,  // 6c_1 good 90
  0.192409082031,  // 6c_2 good 90
  0.184275976562,  // 7c_1 good 90
  0.178534960937,  // 8c_1 good 90
  0.178534960937,  // 8c_2 good 90
  0.175664453125,  // 8c_3 good 90
  0.175664453125,  // 8c_4 good 90
  0.174069726562}; // 8c_5 good 90

double zvalue_95 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 95
  0.207529083252,  // 3_1 good 95
  0.203492431641,  // 4_1 good 95
  0.200781396484,  // 5_1 good 95
  0.199346142578,  // 5_2 good 95
  0.194482226562,  // 6_1 good 95
  0.193924072266,  // 6_2 good 95
  0.193365917969,  // 6_3 good 95
  0.191292773438,  // 7_1 good 95
  0.188980419922,  // 7_2 good 95
  0.188741210938,  // 7_3 good 95
  0.188262792969,  // 7_4 good 95
  0.186827539062,  // 7_5 good 95
  0.186827539062,  // 7_6 good 95
  0.186827539062,  // 7_7 good 95
  0.182202832031,  // 8_1 good 95
  0.181724414062,  // 8_2 good 95
  0.181724414062,  // 8_3 good 95
  0.181485205078,  // 8_4 good 95
  0.180608105469,  // 8_5 good 95
  0.179651269531,  // 8_6 good 95
  0.180608105469,  // 8_7 good 95
  0.179172851562,  // 8_8 good 95
  0.180289160156,  // 8_9 good 95
  0.179172851562,  // 8_10 good 95
  0.179172851562,  // 8_11 good 95
  0.178056542969,  // 8_12 good 95
  0.179651269531,  // 8_13 good 95
  0.178534960937,  // 8_14 good 95
  0.177897070312,  // 8_15 good 95
  0.180289160156,  // 8_16 good 95
  0.179332324219,  // 8_17 good 95
  0.179651269531,  // 8_18 good 95
  0.194083544922,  // 8_19 good 95
  0.190814355469,  // 8_20 good 95
  0.188262792969,  // 8_21 good 95
  UNKNOWN_Z,  // 9_1 unknown 95
  UNKNOWN_Z,  // 9_2 unknown 95
  UNKNOWN_Z,  // 9_3 unknown 95
  UNKNOWN_Z,  // 9_4 unknown 95
  UNKNOWN_Z,  // 9_5 unknown 95
  UNKNOWN_Z,  // 9_6 unknown 95
  UNKNOWN_Z,  // 9_7 unknown 95
  UNKNOWN_Z,  // 9_8 unknown 95
  UNKNOWN_Z,  // 9_9 unknown 95
  UNKNOWN_Z,  // 9_10 unknown 95
  UNKNOWN_Z,  // 9_11 unknown 95
  UNKNOWN_Z,  // 9_12 unknown 95
  UNKNOWN_Z,  // 9_13 unknown 95
  UNKNOWN_Z,  // 9_14 unknown 95
  UNKNOWN_Z,  // 9_15 unknown 95
  UNKNOWN_Z,  // 9_16 unknown 95
  UNKNOWN_Z,  // 9_17 unknown 95
  UNKNOWN_Z,  // 9_18 unknown 95
  UNKNOWN_Z,  // 9_19 unknown 95
  UNKNOWN_Z,  // 9_20 unknown 95
  UNKNOWN_Z,  // 9_21 unknown 95
  UNKNOWN_Z,  // 9_22 unknown 95
  UNKNOWN_Z,  // 9_23 unknown 95
  UNKNOWN_Z,  // 9_24 unknown 95
  UNKNOWN_Z,  // 9_25 unknown 95
  UNKNOWN_Z,  // 9_26 unknown 95
  UNKNOWN_Z,  // 9_27 unknown 95
  UNKNOWN_Z,  // 9_28 unknown 95
  UNKNOWN_Z,  // 9_29 unknown 95
  UNKNOWN_Z,  // 9_30 unknown 95
  UNKNOWN_Z,  // 9_31 unknown 95
  UNKNOWN_Z,  // 9_32 unknown 95
  UNKNOWN_Z,  // 9_33 unknown 95
  UNKNOWN_Z,  // 9_34 unknown 95
  UNKNOWN_Z,  // 9_35 unknown 95
  UNKNOWN_Z,  // 9_36 unknown 95
  UNKNOWN_Z,  // 9_37 unknown 95
  UNKNOWN_Z,  // 9_38 unknown 95
  UNKNOWN_Z,  // 9_39 unknown 95
  UNKNOWN_Z,  // 9_40 unknown 95
  UNKNOWN_Z,  // 9_41 unknown 95
  UNKNOWN_Z,  // 9_42 unknown 95
  UNKNOWN_Z,  // 9_43 unknown 95
  UNKNOWN_Z,  // 9_44 unknown 95
  UNKNOWN_Z,  // 9_45 unknown 95
  UNKNOWN_Z,  // 9_46 unknown 95
  UNKNOWN_Z,  // 9_47 unknown 95
  UNKNOWN_Z,  // 9_48 unknown 95
  UNKNOWN_Z,  // 9_49 unknown 95
  0.194402490234,  // 6c_1 good 95
  0.194402490234,  // 6c_2 good 95
  0.187305957031,  // 7c_1 good 95
  0.182202832031,  // 8c_1 good 95
  0.182202832031,  // 8c_2 good 95
  0.179970214844,  // 8c_3 good 95
  0.179970214844,  // 8c_4 good 95
  0.178534960937}; // 8c_5 good 95

double zvalue_100 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 100
  0.20799753418,  // 3_1 good 100
  0.204489135742,  // 4_1 good 100
  0.201897705078,  // 5_1 good 100
  0.200735282898,  // 5_2 good 100
  0.196433792114,  // 6_1 good 100
  0.195954589844,  // 6_2 good 100
  0.195279846191,  // 6_3 good 100
  0.193551025391,  // 7_1 good 100
  0.19146463871,  // 7_2 good 100
  0.191174316406,  // 7_3 good 100
  0.190758056641,  // 7_4 good 100
  0.1898097229,  // 7_5 good 100
  0.189415283203,  // 7_6 good 100
  0.189348144531,  // 7_7 good 100
  0.185668945312,  // 8_1 good 100
  0.185104980469,  // 8_2 good 100
  0.185240097046,  // 8_3 good 100
  0.184903564453,  // 8_4 good 100
  0.184044189453,  // 8_5 good 100
  0.1833543396,  // 8_6 good 100
  0.183977050781,  // 8_7 good 100
  0.183090820313,  // 8_8 good 100
  0.183795776367,  // 8_9 good 100
  0.182902832031,  // 8_10 good 100
  0.183063964844,  // 8_11 good 100
  0.181979675293,  // 8_12 good 100
  0.183201599121,  // 8_13 good 100
  0.182124023438,  // 8_14 good 100
  0.182026672363,  // 8_15 good 100
  0.183674926758,  // 8_16 good 100
  0.183010253906,  // 8_17 good 100
  0.183255310059,  // 8_18 good 100
  0.195981445313,  // 8_19 good 100
  0.193243026733,  // 8_20 good 100
  0.190959472656,  // 8_21 good 100
  0.181724414062,  // 9_1 good 100
  0.179172851562,  // 9_2 good 100
  0.178534960937,  // 9_3 good 100
  0.178534960937,  // 9_4 good 100
  0.177897070312,  // 9_5 good 100
  0.176621289062,  // 9_6 good 100
  0.175983398437,  // 9_7 good 100
  0.175983398437,  // 9_8 good 100
  0.175983398437,  // 9_9 good 100
  0.175983398437,  // 9_10 good 100
  0.175983398437,  // 9_11 good 100
  0.175983398437,  // 9_12 good 100
  0.175664453125,  // 9_13 good 100
  0.176621289062,  // 9_14 good 100
  0.174069726562,  // 9_15 good 100
  0.175186035156,  // 9_16 good 100
  0.175664453125,  // 9_17 good 100
  0.174069726562,  // 9_18 good 100
  0.174867089844,  // 9_19 good 100
  0.175504980469,  // 9_20 good 100
  0.174548144531,  // 9_21 good 100
  0.174867089844,  // 9_22 good 100
  0.174069726562,  // 9_23 good 100
  0.174388671875,  // 9_24 good 100
  0.173431835938,  // 9_25 good 100
  0.175186035156,  // 9_26 good 100
  0.174388671875,  // 9_27 good 100
  0.174069726562,  // 9_28 good 100
  0.175504980469,  // 9_29 good 100
  0.174069726562,  // 9_30 good 100
  0.174388671875,  // 9_31 good 100
  0.175664453125,  // 9_32 good 100
  0.175504980469,  // 9_33 good 100
  0.175504980469,  // 9_34 good 100
  0.177099707031,  // 9_35 good 100
  0.174867089844,  // 9_36 good 100
  0.174867089844,  // 9_37 good 100
  0.174069726562,  // 9_38 good 100
  0.174867089844,  // 9_39 good 100
  0.180608105469,  // 9_40 good 100
  0.175983398437,  // 9_41 good 100
  0.188741210938,  // 9_42 good 100
  0.186827539062,  // 9_43 good 100
  0.185711230469,  // 9_44 good 100
  0.183877294922,  // 9_45 good 100
  0.189379101562,  // 9_46 good 100
  0.185711230469,  // 9_47 good 100
  0.183638085938,  // 9_48 good 100
  0.184275976562,  // 9_49 good 100
  0.196146459961,  // 6c_1 good 100
  0.196206298828,  // 6c_2 good 100
  0.189943164063,  // 7c_1 good 100
  0.185395410156,  // 8c_1 good 100
  0.185554980469,  // 8c_2 good 100
  0.183560351562,  // 8c_3 good 100
  0.183560351562,  // 8c_4 good 100
  0.182283789063}; // 8c_5 good 100

double zvalue_105 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 105
  0.20849588623,  // 3_1 good 105
  0.205246630859,  // 4_1 good 105
  0.203014013672,  // 5_1 good 105
  0.201857836914,  // 5_2 good 105
  0.198070361328,  // 6_1 good 105
  0.197512207031,  // 6_2 good 105
  0.196954052734,  // 6_3 good 105
  0.195518798828,  // 7_1 good 105
  0.193684863281,  // 7_2 good 105
  0.193365917969,  // 7_3 good 105
  0.192967236328,  // 7_4 good 105
  0.192090136719,  // 7_5 good 105
  0.191691455078,  // 7_6 good 105
  0.191691455078,  // 7_7 good 105
  0.188581738281,  // 8_1 good 105
  0.187943847656,  // 8_2 good 105
  0.187943847656,  // 8_3 good 105
  0.187704638672,  // 8_4 good 105
  0.186827539062,  // 8_5 good 105
  0.186189648438,  // 8_6 good 105
  0.186827539062,  // 8_7 good 105
  0.186189648438,  // 8_8 good 105
  0.186827539062,  // 8_9 good 105
  0.186030175781,  // 8_10 good 105
  0.186189648438,  // 8_11 good 105
  0.185312548828,  // 8_12 good 105
  0.186189648438,  // 8_13 good 105
  0.185392285156,  // 8_14 good 105
  0.185073339844,  // 8_15 good 105
  0.186588330078,  // 8_16 good 105
  0.186030175781,  // 8_17 good 105
  0.186189648438,  // 8_18 good 105
  0.197831152344,  // 8_19 good 105
  0.195199853516,  // 8_20 good 105
  0.193046972656,  // 8_21 good 105
  UNKNOWN_Z,  // 9_1 unknown 105
  UNKNOWN_Z,  // 9_2 unknown 105
  UNKNOWN_Z,  // 9_3 unknown 105
  UNKNOWN_Z,  // 9_4 unknown 105
  UNKNOWN_Z,  // 9_5 unknown 105
  UNKNOWN_Z,  // 9_6 unknown 105
  UNKNOWN_Z,  // 9_7 unknown 105
  UNKNOWN_Z,  // 9_8 unknown 105
  UNKNOWN_Z,  // 9_9 unknown 105
  UNKNOWN_Z,  // 9_10 unknown 105
  UNKNOWN_Z,  // 9_11 unknown 105
  UNKNOWN_Z,  // 9_12 unknown 105
  UNKNOWN_Z,  // 9_13 unknown 105
  UNKNOWN_Z,  // 9_14 unknown 105
  UNKNOWN_Z,  // 9_15 unknown 105
  UNKNOWN_Z,  // 9_16 unknown 105
  UNKNOWN_Z,  // 9_17 unknown 105
  UNKNOWN_Z,  // 9_18 unknown 105
  UNKNOWN_Z,  // 9_19 unknown 105
  UNKNOWN_Z,  // 9_20 unknown 105
  UNKNOWN_Z,  // 9_21 unknown 105
  UNKNOWN_Z,  // 9_22 unknown 105
  UNKNOWN_Z,  // 9_23 unknown 105
  UNKNOWN_Z,  // 9_24 unknown 105
  UNKNOWN_Z,  // 9_25 unknown 105
  UNKNOWN_Z,  // 9_26 unknown 105
  UNKNOWN_Z,  // 9_27 unknown 105
  UNKNOWN_Z,  // 9_28 unknown 105
  UNKNOWN_Z,  // 9_29 unknown 105
  UNKNOWN_Z,  // 9_30 unknown 105
  UNKNOWN_Z,  // 9_31 unknown 105
  UNKNOWN_Z,  // 9_32 unknown 105
  UNKNOWN_Z,  // 9_33 unknown 105
  UNKNOWN_Z,  // 9_34 unknown 105
  UNKNOWN_Z,  // 9_35 unknown 105
  UNKNOWN_Z,  // 9_36 unknown 105
  UNKNOWN_Z,  // 9_37 unknown 105
  UNKNOWN_Z,  // 9_38 unknown 105
  UNKNOWN_Z,  // 9_39 unknown 105
  UNKNOWN_Z,  // 9_40 unknown 105
  UNKNOWN_Z,  // 9_41 unknown 105
  UNKNOWN_Z,  // 9_42 unknown 105
  UNKNOWN_Z,  // 9_43 unknown 105
  UNKNOWN_Z,  // 9_44 unknown 105
  UNKNOWN_Z,  // 9_45 unknown 105
  UNKNOWN_Z,  // 9_46 unknown 105
  UNKNOWN_Z,  // 9_47 unknown 105
  UNKNOWN_Z,  // 9_48 unknown 105
  UNKNOWN_Z,  // 9_49 unknown 105
  0.197910888672,  // 6c_1 good 105
  0.197751416016,  // 6c_2 good 105
  0.192090136719,  // 7c_1 good 105
  0.188262792969,  // 8c_1 good 105
  0.188262792969,  // 8c_2 good 105
  0.186588330078,  // 8c_3 good 105
  0.186588330078,  // 8c_4 good 105
  0.185392285156}; // 8c_5 good 105

double zvalue_110 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 110
  0.208894567871,  // 3_1 good 110
  0.20594432373,  // 4_1 good 110
  0.203811376953,  // 5_1 good 110
  0.202854541016,  // 5_2 good 110
  0.199505615234,  // 6_1 good 110
  0.199027197266,  // 6_2 good 110
  0.198469042969,  // 6_3 good 110
  0.197193261719,  // 7_1 good 110
  0.195518798828,  // 7_2 good 110
  0.195199853516,  // 7_3 good 110
  0.194960644531,  // 7_4 good 110
  0.194083544922,  // 7_5 good 110
  0.193684863281,  // 7_6 good 110
  0.193684863281,  // 7_7 good 110
  0.190814355469,  // 8_1 good 110
  0.190415673828,  // 8_2 good 110
  0.190415673828,  // 8_3 good 110
  0.190096728516,  // 8_4 good 110
  0.189379101562,  // 8_5 good 110
  0.188980419922,  // 8_6 good 110
  0.189379101562,  // 8_7 good 110
  0.188741210938,  // 8_8 good 110
  0.189299365234,  // 8_9 good 110
  0.188581738281,  // 8_10 good 110
  0.188741210938,  // 8_11 good 110
  0.187864111328,  // 8_12 good 110
  0.188741210938,  // 8_13 good 110
  0.187943847656,  // 8_14 good 110
  0.187864111328,  // 8_15 good 110
  0.188980419922,  // 8_16 good 110
  0.188581738281,  // 8_17 good 110
  0.188581738281,  // 8_18 good 110
  0.199027197266,  // 8_19 good 110
  0.196794580078,  // 8_20 good 110
  0.195279589844,  // 8_21 good 110
  UNKNOWN_Z,  // 9_1 unknown 110
  UNKNOWN_Z,  // 9_2 unknown 110
  UNKNOWN_Z,  // 9_3 unknown 110
  UNKNOWN_Z,  // 9_4 unknown 110
  UNKNOWN_Z,  // 9_5 unknown 110
  UNKNOWN_Z,  // 9_6 unknown 110
  UNKNOWN_Z,  // 9_7 unknown 110
  UNKNOWN_Z,  // 9_8 unknown 110
  UNKNOWN_Z,  // 9_9 unknown 110
  UNKNOWN_Z,  // 9_10 unknown 110
  UNKNOWN_Z,  // 9_11 unknown 110
  UNKNOWN_Z,  // 9_12 unknown 110
  UNKNOWN_Z,  // 9_13 unknown 110
  UNKNOWN_Z,  // 9_14 unknown 110
  UNKNOWN_Z,  // 9_15 unknown 110
  UNKNOWN_Z,  // 9_16 unknown 110
  UNKNOWN_Z,  // 9_17 unknown 110
  UNKNOWN_Z,  // 9_18 unknown 110
  UNKNOWN_Z,  // 9_19 unknown 110
  UNKNOWN_Z,  // 9_20 unknown 110
  UNKNOWN_Z,  // 9_21 unknown 110
  UNKNOWN_Z,  // 9_22 unknown 110
  UNKNOWN_Z,  // 9_23 unknown 110
  UNKNOWN_Z,  // 9_24 unknown 110
  UNKNOWN_Z,  // 9_25 unknown 110
  UNKNOWN_Z,  // 9_26 unknown 110
  UNKNOWN_Z,  // 9_27 unknown 110
  UNKNOWN_Z,  // 9_28 unknown 110
  UNKNOWN_Z,  // 9_29 unknown 110
  UNKNOWN_Z,  // 9_30 unknown 110
  UNKNOWN_Z,  // 9_31 unknown 110
  UNKNOWN_Z,  // 9_32 unknown 110
  UNKNOWN_Z,  // 9_33 unknown 110
  UNKNOWN_Z,  // 9_34 unknown 110
  UNKNOWN_Z,  // 9_35 unknown 110
  UNKNOWN_Z,  // 9_36 unknown 110
  UNKNOWN_Z,  // 9_37 unknown 110
  UNKNOWN_Z,  // 9_38 unknown 110
  UNKNOWN_Z,  // 9_39 unknown 110
  UNKNOWN_Z,  // 9_40 unknown 110
  UNKNOWN_Z,  // 9_41 unknown 110
  UNKNOWN_Z,  // 9_42 unknown 110
  UNKNOWN_Z,  // 9_43 unknown 110
  UNKNOWN_Z,  // 9_44 unknown 110
  UNKNOWN_Z,  // 9_45 unknown 110
  UNKNOWN_Z,  // 9_46 unknown 110
  UNKNOWN_Z,  // 9_47 unknown 110
  UNKNOWN_Z,  // 9_48 unknown 110
  UNKNOWN_Z,  // 9_49 unknown 110
  0.199186669922,  // 6c_1 good 110
  0.199186669922,  // 6c_2 good 110
  0.194083544922,  // 7c_1 good 110
  0.190495410156,  // 8c_1 good 110
  0.190495410156,  // 8c_2 good 110
  0.189139892578,  // 8c_3 good 110
  0.189139892578,  // 8c_4 good 110
  0.188262792969}; // 8c_5 good 110

double zvalue_115 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 115
  0.209243414307,  // 3_1 good 115
  0.206612115479,  // 4_1 good 115
  0.204608740234,  // 5_1 good 115
  0.203851245117,  // 5_2 good 115
  0.200661791992,  // 6_1 good 115
  0.200302978516,  // 6_2 good 115
  0.199744824219,  // 6_3 good 115
  0.198469042969,  // 7_1 good 115
  0.197113525391,  // 7_2 good 115
  0.196794580078,  // 7_3 good 115
  0.196475634766,  // 7_4 good 115
  0.195917480469,  // 7_5 good 115
  0.195518798828,  // 7_6 good 115
  0.195359326172,  // 7_7 good 115
  0.192807763672,  // 8_1 good 115
  0.192409082031,  // 8_2 good 115
  0.192409082031,  // 8_3 good 115
  0.192090136719,  // 8_4 good 115
  0.191691455078,  // 8_5 good 115
  0.191133300781,  // 8_6 good 115
  0.191531982422,  // 8_7 good 115
  0.191133300781,  // 8_8 good 115
  0.191292773438,  // 8_9 good 115
  0.190814355469,  // 8_10 good 115
  0.190814355469,  // 8_11 good 115
  0.190256201172,  // 8_12 good 115
  0.190814355469,  // 8_13 good 115
  0.190256201172,  // 8_14 good 115
  0.190176464844,  // 8_15 good 115
  0.191133300781,  // 8_16 good 115
  0.190495410156,  // 8_17 good 115
  0.190495410156,  // 8_18 good 115
  0.200302978516,  // 8_19 good 115
  0.198389306641,  // 8_20 good 115
  0.196834448242,  // 8_21 good 115
  UNKNOWN_Z,  // 9_1 unknown 115
  UNKNOWN_Z,  // 9_2 unknown 115
  UNKNOWN_Z,  // 9_3 unknown 115
  UNKNOWN_Z,  // 9_4 unknown 115
  UNKNOWN_Z,  // 9_5 unknown 115
  UNKNOWN_Z,  // 9_6 unknown 115
  UNKNOWN_Z,  // 9_7 unknown 115
  UNKNOWN_Z,  // 9_8 unknown 115
  UNKNOWN_Z,  // 9_9 unknown 115
  UNKNOWN_Z,  // 9_10 unknown 115
  UNKNOWN_Z,  // 9_11 unknown 115
  UNKNOWN_Z,  // 9_12 unknown 115
  UNKNOWN_Z,  // 9_13 unknown 115
  UNKNOWN_Z,  // 9_14 unknown 115
  UNKNOWN_Z,  // 9_15 unknown 115
  UNKNOWN_Z,  // 9_16 unknown 115
  UNKNOWN_Z,  // 9_17 unknown 115
  UNKNOWN_Z,  // 9_18 unknown 115
  UNKNOWN_Z,  // 9_19 unknown 115
  UNKNOWN_Z,  // 9_20 unknown 115
  UNKNOWN_Z,  // 9_21 unknown 115
  UNKNOWN_Z,  // 9_22 unknown 115
  UNKNOWN_Z,  // 9_23 unknown 115
  UNKNOWN_Z,  // 9_24 unknown 115
  UNKNOWN_Z,  // 9_25 unknown 115
  UNKNOWN_Z,  // 9_26 unknown 115
  UNKNOWN_Z,  // 9_27 unknown 115
  UNKNOWN_Z,  // 9_28 unknown 115
  UNKNOWN_Z,  // 9_29 unknown 115
  UNKNOWN_Z,  // 9_30 unknown 115
  UNKNOWN_Z,  // 9_31 unknown 115
  UNKNOWN_Z,  // 9_32 unknown 115
  UNKNOWN_Z,  // 9_33 unknown 115
  UNKNOWN_Z,  // 9_34 unknown 115
  UNKNOWN_Z,  // 9_35 unknown 115
  UNKNOWN_Z,  // 9_36 unknown 115
  UNKNOWN_Z,  // 9_37 unknown 115
  UNKNOWN_Z,  // 9_38 unknown 115
  UNKNOWN_Z,  // 9_39 unknown 115
  UNKNOWN_Z,  // 9_40 unknown 115
  UNKNOWN_Z,  // 9_41 unknown 115
  UNKNOWN_Z,  // 9_42 unknown 115
  UNKNOWN_Z,  // 9_43 unknown 115
  UNKNOWN_Z,  // 9_44 unknown 115
  UNKNOWN_Z,  // 9_45 unknown 115
  UNKNOWN_Z,  // 9_46 unknown 115
  UNKNOWN_Z,  // 9_47 unknown 115
  UNKNOWN_Z,  // 9_48 unknown 115
  UNKNOWN_Z,  // 9_49 unknown 115
  0.200302978516,  // 6c_1 good 115
  0.200302978516,  // 6c_2 good 115
  0.195598535156,  // 7c_1 good 115
  0.192728027344,  // 8c_1 good 115
  0.192728027344,  // 8c_2 good 115
  0.191292773438,  // 8c_3 good 115
  0.191292773438,  // 8c_4 good 115
  0.190415673828}; // 8c_5 good 115

double zvalue_120 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 120
  0.209612194824,  // 3_1 good 120
  0.207078074646,  // 4_1 good 120
  0.20534630127,  // 5_1 good 120
  0.20456887207,  // 5_2 good 120
  0.201578759766,  // 6_1 good 120
  0.201299682617,  // 6_2 good 120
  0.200940869141,  // 6_3 good 120
  0.199744824219,  // 7_1 good 120
  0.198389306641,  // 7_2 good 120
  0.198229833984,  // 7_3 good 120
  0.197751416016,  // 7_4 good 120
  0.197193261719,  // 7_5 good 120
  0.196954052734,  // 7_6 good 120
  0.196794580078,  // 7_7 good 120
  0.194641699219,  // 8_1 good 120
  0.194243017578,  // 8_2 good 120
  0.194243017578,  // 8_3 good 120
  0.193924072266,  // 8_4 good 120
  0.193365917969,  // 8_5 good 120
  0.193126708984,  // 8_6 good 120
  0.193365917969,  // 8_7 good 120
  0.192967236328,  // 8_8 good 120
  0.193365917969,  // 8_9 good 120
  0.192648291016,  // 8_10 good 120
  0.192807763672,  // 8_11 good 120
  0.192409082031,  // 8_12 good 120
  0.192648291016,  // 8_13 good 120
  0.192409082031,  // 8_14 good 120
  0.192090136719,  // 8_15 good 120
  0.192927368164,  // 8_16 good 120
  0.192409082031,  // 8_17 good 120
  0.192409082031,  // 8_18 good 120
  0.201419287109,  // 8_19 good 120
  0.199505615234,  // 8_20 good 120
  0.198229833984,  // 8_21 good 120
  UNKNOWN_Z,  // 9_1 unknown 120
  UNKNOWN_Z,  // 9_2 unknown 120
  UNKNOWN_Z,  // 9_3 unknown 120
  UNKNOWN_Z,  // 9_4 unknown 120
  UNKNOWN_Z,  // 9_5 unknown 120
  UNKNOWN_Z,  // 9_6 unknown 120
  UNKNOWN_Z,  // 9_7 unknown 120
  UNKNOWN_Z,  // 9_8 unknown 120
  UNKNOWN_Z,  // 9_9 unknown 120
  UNKNOWN_Z,  // 9_10 unknown 120
  UNKNOWN_Z,  // 9_11 unknown 120
  UNKNOWN_Z,  // 9_12 unknown 120
  UNKNOWN_Z,  // 9_13 unknown 120
  UNKNOWN_Z,  // 9_14 unknown 120
  UNKNOWN_Z,  // 9_15 unknown 120
  UNKNOWN_Z,  // 9_16 unknown 120
  UNKNOWN_Z,  // 9_17 unknown 120
  UNKNOWN_Z,  // 9_18 unknown 120
  UNKNOWN_Z,  // 9_19 unknown 120
  UNKNOWN_Z,  // 9_20 unknown 120
  UNKNOWN_Z,  // 9_21 unknown 120
  UNKNOWN_Z,  // 9_22 unknown 120
  UNKNOWN_Z,  // 9_23 unknown 120
  UNKNOWN_Z,  // 9_24 unknown 120
  UNKNOWN_Z,  // 9_25 unknown 120
  UNKNOWN_Z,  // 9_26 unknown 120
  UNKNOWN_Z,  // 9_27 unknown 120
  UNKNOWN_Z,  // 9_28 unknown 120
  UNKNOWN_Z,  // 9_29 unknown 120
  UNKNOWN_Z,  // 9_30 unknown 120
  UNKNOWN_Z,  // 9_31 unknown 120
  UNKNOWN_Z,  // 9_32 unknown 120
  UNKNOWN_Z,  // 9_33 unknown 120
  UNKNOWN_Z,  // 9_34 unknown 120
  UNKNOWN_Z,  // 9_35 unknown 120
  UNKNOWN_Z,  // 9_36 unknown 120
  UNKNOWN_Z,  // 9_37 unknown 120
  UNKNOWN_Z,  // 9_38 unknown 120
  UNKNOWN_Z,  // 9_39 unknown 120
  UNKNOWN_Z,  // 9_40 unknown 120
  UNKNOWN_Z,  // 9_41 unknown 120
  UNKNOWN_Z,  // 9_42 unknown 120
  UNKNOWN_Z,  // 9_43 unknown 120
  UNKNOWN_Z,  // 9_44 unknown 120
  UNKNOWN_Z,  // 9_45 unknown 120
  UNKNOWN_Z,  // 9_46 unknown 120
  UNKNOWN_Z,  // 9_47 unknown 120
  UNKNOWN_Z,  // 9_48 unknown 120
  UNKNOWN_Z,  // 9_49 unknown 120
  0.201339550781,  // 6c_1 good 120
  0.201339550781,  // 6c_2 good 120
  0.197193261719,  // 7c_1 good 120
  0.194402490234,  // 8c_1 good 120
  0.194402490234,  // 8c_2 good 120
  0.193126708984,  // 8c_3 good 120
  0.193126708984,  // 8c_4 good 120
  0.192409082031}; // 8c_5 good 120

double zvalue_124 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 124
  0.209758390033,  // 3_1 good 124
  0.207480938721,  // 4_1 good 124
  0.205780517578,  // 5_1 good 124
  0.205122290039,  // 5_2 good 124
  0.202349755859,  // 6_1 good 124
  0.202030615234,  // 6_2 good 124
  0.201651635742,  // 6_3 good 124
  0.200594482422,  // 7_1 good 124
  0.199317919922,  // 7_2 good 124
  0.199118457031,  // 7_3 good 124
  0.198839208984,  // 7_4 good 124
  0.198320605469,  // 7_5 good 124
  0.198041357422,  // 7_6 good 124
  0.197841894531,  // 7_7 good 124
  0.195927050781,  // 8_1 good 124
  0.195448339844,  // 8_2 good 124
  0.195568017578,  // 8_3 good 124
  0.195288769531,  // 8_4 good 124
  0.194730273438,  // 8_5 good 124
  0.194490917969,  // 8_6 good 124
  0.194730273438,  // 8_7 good 124
  0.194211669922,  // 8_8 good 124
  0.194490917969,  // 8_9 good 124
  0.194012207031,  // 8_10 good 124
  0.194171777344,  // 8_11 good 124
  0.193653173828,  // 8_12 good 124
  0.194012207031,  // 8_13 good 124
  0.193653173828,  // 8_14 good 124
  0.193453710937,  // 8_15 good 124
  0.194131884766,  // 8_16 good 124
  0.193852636719,  // 8_17 good 124
  0.193772851562,  // 8_18 good 124
  0.202150292969,  // 8_19 good 124
  0.200395019531,  // 8_20 good 124
  0.199118457031,  // 8_21 good 124
  UNKNOWN_Z,  // 9_1 unknown 124
  UNKNOWN_Z,  // 9_2 unknown 124
  UNKNOWN_Z,  // 9_3 unknown 124
  UNKNOWN_Z,  // 9_4 unknown 124
  UNKNOWN_Z,  // 9_5 unknown 124
  UNKNOWN_Z,  // 9_6 unknown 124
  UNKNOWN_Z,  // 9_7 unknown 124
  UNKNOWN_Z,  // 9_8 unknown 124
  UNKNOWN_Z,  // 9_9 unknown 124
  UNKNOWN_Z,  // 9_10 unknown 124
  UNKNOWN_Z,  // 9_11 unknown 124
  UNKNOWN_Z,  // 9_12 unknown 124
  UNKNOWN_Z,  // 9_13 unknown 124
  UNKNOWN_Z,  // 9_14 unknown 124
  UNKNOWN_Z,  // 9_15 unknown 124
  UNKNOWN_Z,  // 9_16 unknown 124
  UNKNOWN_Z,  // 9_17 unknown 124
  UNKNOWN_Z,  // 9_18 unknown 124
  UNKNOWN_Z,  // 9_19 unknown 124
  UNKNOWN_Z,  // 9_20 unknown 124
  UNKNOWN_Z,  // 9_21 unknown 124
  UNKNOWN_Z,  // 9_22 unknown 124
  UNKNOWN_Z,  // 9_23 unknown 124
  UNKNOWN_Z,  // 9_24 unknown 124
  UNKNOWN_Z,  // 9_25 unknown 124
  UNKNOWN_Z,  // 9_26 unknown 124
  UNKNOWN_Z,  // 9_27 unknown 124
  UNKNOWN_Z,  // 9_28 unknown 124
  UNKNOWN_Z,  // 9_29 unknown 124
  UNKNOWN_Z,  // 9_30 unknown 124
  UNKNOWN_Z,  // 9_31 unknown 124
  UNKNOWN_Z,  // 9_32 unknown 124
  UNKNOWN_Z,  // 9_33 unknown 124
  UNKNOWN_Z,  // 9_34 unknown 124
  UNKNOWN_Z,  // 9_35 unknown 124
  UNKNOWN_Z,  // 9_36 unknown 124
  UNKNOWN_Z,  // 9_37 unknown 124
  UNKNOWN_Z,  // 9_38 unknown 124
  UNKNOWN_Z,  // 9_39 unknown 124
  UNKNOWN_Z,  // 9_40 unknown 124
  UNKNOWN_Z,  // 9_41 unknown 124
  UNKNOWN_Z,  // 9_42 unknown 124
  UNKNOWN_Z,  // 9_43 unknown 124
  UNKNOWN_Z,  // 9_44 unknown 124
  UNKNOWN_Z,  // 9_45 unknown 124
  UNKNOWN_Z,  // 9_46 unknown 124
  UNKNOWN_Z,  // 9_47 unknown 124
  UNKNOWN_Z,  // 9_48 unknown 124
  UNKNOWN_Z,  // 9_49 unknown 124
  0.202030615234,  // 6c_1 good 124
  0.201990722656,  // 6c_2 good 124
  0.198229833984,  // 7c_1 good 124
  0.195518798828,  // 8c_1 good 124
  0.195518798828,  // 8c_2 good 124
  0.194402490234,  // 8c_3 good 124
  0.194402490234,  // 8c_4 good 124
  0.193684863281}; // 8c_5 good 124

double zvalue_125 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 125
  0.209784126282,  // 3_1 good 125
  0.207628753662,  // 4_1 good 125
  0.20594432373,  // 5_1 good 125
  0.205207459331,  // 5_2 good 125
  0.202575463867,  // 6_1 good 125
  0.202208557129,  // 6_2 good 125
  0.201793766022,  // 6_3 good 125
  0.200781396484,  // 7_1 good 125
  0.199665087891,  // 7_2 good 125
  0.199386010742,  // 7_3 good 125
  0.199027197266,  // 7_4 good 125
  0.198542785645,  // 7_5 good 125
  0.198229833984,  // 7_6 good 125
  0.198229833984,  // 7_7 good 125
  0.196236425781,  // 8_1 good 125
  0.195679321289,  // 8_2 good 125
  0.195818371773,  // 8_3 good 125
  0.195549158156,  // 8_4 good 125
  0.195061645508,  // 8_5 good 125
  0.194641699219,  // 8_6 good 125
  0.194959259033,  // 8_7 good 125
  0.194522015154,  // 8_8 good 125
  0.194960644531,  // 8_9 good 125
  0.194356584549,  // 8_10 good 125
  0.194464111328,  // 8_11 good 125
  0.193967285156,  // 8_12 good 125
  0.194403686523,  // 8_13 good 125
  0.193958892822,  // 8_14 good 125
  0.193850631714,  // 8_15 good 125
  0.194402490234,  // 8_16 good 125
  0.19418926239,  // 8_17 good 125
  0.194178771973,  // 8_18 good 125
  0.202296386719,  // 8_19 good 125
  0.200582055664,  // 8_20 good 125
  0.199346142578,  // 8_21 good 125
  UNKNOWN_Z,  // 9_1 unknown 125
  UNKNOWN_Z,  // 9_2 unknown 125
  UNKNOWN_Z,  // 9_3 unknown 125
  UNKNOWN_Z,  // 9_4 unknown 125
  UNKNOWN_Z,  // 9_5 unknown 125
  UNKNOWN_Z,  // 9_6 unknown 125
  UNKNOWN_Z,  // 9_7 unknown 125
  UNKNOWN_Z,  // 9_8 unknown 125
  UNKNOWN_Z,  // 9_9 unknown 125
  UNKNOWN_Z,  // 9_10 unknown 125
  UNKNOWN_Z,  // 9_11 unknown 125
  UNKNOWN_Z,  // 9_12 unknown 125
  UNKNOWN_Z,  // 9_13 unknown 125
  UNKNOWN_Z,  // 9_14 unknown 125
  UNKNOWN_Z,  // 9_15 unknown 125
  UNKNOWN_Z,  // 9_16 unknown 125
  UNKNOWN_Z,  // 9_17 unknown 125
  UNKNOWN_Z,  // 9_18 unknown 125
  UNKNOWN_Z,  // 9_19 unknown 125
  UNKNOWN_Z,  // 9_20 unknown 125
  UNKNOWN_Z,  // 9_21 unknown 125
  UNKNOWN_Z,  // 9_22 unknown 125
  UNKNOWN_Z,  // 9_23 unknown 125
  UNKNOWN_Z,  // 9_24 unknown 125
  UNKNOWN_Z,  // 9_25 unknown 125
  UNKNOWN_Z,  // 9_26 unknown 125
  UNKNOWN_Z,  // 9_27 unknown 125
  UNKNOWN_Z,  // 9_28 unknown 125
  UNKNOWN_Z,  // 9_29 unknown 125
  UNKNOWN_Z,  // 9_30 unknown 125
  UNKNOWN_Z,  // 9_31 unknown 125
  UNKNOWN_Z,  // 9_32 unknown 125
  UNKNOWN_Z,  // 9_33 unknown 125
  UNKNOWN_Z,  // 9_34 unknown 125
  UNKNOWN_Z,  // 9_35 unknown 125
  UNKNOWN_Z,  // 9_36 unknown 125
  UNKNOWN_Z,  // 9_37 unknown 125
  UNKNOWN_Z,  // 9_38 unknown 125
  UNKNOWN_Z,  // 9_39 unknown 125
  UNKNOWN_Z,  // 9_40 unknown 125
  UNKNOWN_Z,  // 9_41 unknown 125
  UNKNOWN_Z,  // 9_42 unknown 125
  UNKNOWN_Z,  // 9_43 unknown 125
  UNKNOWN_Z,  // 9_44 unknown 125
  UNKNOWN_Z,  // 9_45 unknown 125
  UNKNOWN_Z,  // 9_46 unknown 125
  UNKNOWN_Z,  // 9_47 unknown 125
  UNKNOWN_Z,  // 9_48 unknown 125
  UNKNOWN_Z,  // 9_49 unknown 125
  0.202216650391,  // 6c_1 good 125
  0.202216650391,  // 6c_2 good 125
  0.198469042969,  // 7c_1 good 125
  0.195837744141,  // 8c_1 good 125
  0.195837744141,  // 8c_2 good 125
  0.194641699219,  // 8c_3 good 125
  0.194641699219,  // 8c_4 good 125
  0.193924072266}; // 8c_5 good 125

double zvalue_130 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 130
  0.21011101408,  // 3_1 good 130
  0.208042385864,  // 4_1 good 130
  0.206522412109,  // 5_1 good 130
  0.205884521484,  // 5_2 good 130
  0.203332958984,  // 6_1 good 130
  0.202974145508,  // 6_2 good 130
  0.202535595703,  // 6_3 good 130
  0.201698364258,  // 7_1 good 130
  0.200621923828,  // 7_2 good 130
  0.200462451172,  // 7_3 good 130
  0.200063769531,  // 7_4 good 130
  0.199744824219,  // 7_5 good 130
  0.199346142578,  // 7_6 good 130
  0.199306274414,  // 7_7 good 130
  0.197512207031,  // 8_1 good 130
  0.197113525391,  // 8_2 good 130
  0.197193261719,  // 8_3 good 130
  0.196954052734,  // 8_4 good 130
  0.196475634766,  // 8_5 good 130
  0.196236425781,  // 8_6 good 130
  0.196475634766,  // 8_7 good 130
  0.195917480469,  // 8_8 good 130
  0.196236425781,  // 8_9 good 130
  0.195917480469,  // 8_10 good 130
  0.195917480469,  // 8_11 good 130
  0.195518798828,  // 8_12 good 130
  0.195917480469,  // 8_13 good 130
  0.195478930664,  // 8_14 good 130
  0.195359326172,  // 8_15 good 130
  0.195917480469,  // 8_16 good 130
  0.195518798828,  // 8_17 good 130
  0.195518798828,  // 8_18 good 130
  0.203133618164,  // 8_19 good 130
  0.201578759766,  // 8_20 good 130
  0.200462451172,  // 8_21 good 130
  UNKNOWN_Z,  // 9_1 unknown 130
  UNKNOWN_Z,  // 9_2 unknown 130
  UNKNOWN_Z,  // 9_3 unknown 130
  UNKNOWN_Z,  // 9_4 unknown 130
  UNKNOWN_Z,  // 9_5 unknown 130
  UNKNOWN_Z,  // 9_6 unknown 130
  UNKNOWN_Z,  // 9_7 unknown 130
  UNKNOWN_Z,  // 9_8 unknown 130
  UNKNOWN_Z,  // 9_9 unknown 130
  UNKNOWN_Z,  // 9_10 unknown 130
  UNKNOWN_Z,  // 9_11 unknown 130
  UNKNOWN_Z,  // 9_12 unknown 130
  UNKNOWN_Z,  // 9_13 unknown 130
  UNKNOWN_Z,  // 9_14 unknown 130
  UNKNOWN_Z,  // 9_15 unknown 130
  UNKNOWN_Z,  // 9_16 unknown 130
  UNKNOWN_Z,  // 9_17 unknown 130
  UNKNOWN_Z,  // 9_18 unknown 130
  UNKNOWN_Z,  // 9_19 unknown 130
  UNKNOWN_Z,  // 9_20 unknown 130
  UNKNOWN_Z,  // 9_21 unknown 130
  UNKNOWN_Z,  // 9_22 unknown 130
  UNKNOWN_Z,  // 9_23 unknown 130
  UNKNOWN_Z,  // 9_24 unknown 130
  UNKNOWN_Z,  // 9_25 unknown 130
  UNKNOWN_Z,  // 9_26 unknown 130
  UNKNOWN_Z,  // 9_27 unknown 130
  UNKNOWN_Z,  // 9_28 unknown 130
  UNKNOWN_Z,  // 9_29 unknown 130
  UNKNOWN_Z,  // 9_30 unknown 130
  UNKNOWN_Z,  // 9_31 unknown 130
  UNKNOWN_Z,  // 9_32 unknown 130
  UNKNOWN_Z,  // 9_33 unknown 130
  UNKNOWN_Z,  // 9_34 unknown 130
  UNKNOWN_Z,  // 9_35 unknown 130
  UNKNOWN_Z,  // 9_36 unknown 130
  UNKNOWN_Z,  // 9_37 unknown 130
  UNKNOWN_Z,  // 9_38 unknown 130
  UNKNOWN_Z,  // 9_39 unknown 130
  UNKNOWN_Z,  // 9_40 unknown 130
  UNKNOWN_Z,  // 9_41 unknown 130
  UNKNOWN_Z,  // 9_42 unknown 130
  UNKNOWN_Z,  // 9_43 unknown 130
  UNKNOWN_Z,  // 9_44 unknown 130
  UNKNOWN_Z,  // 9_45 unknown 130
  UNKNOWN_Z,  // 9_46 unknown 130
  UNKNOWN_Z,  // 9_47 unknown 130
  UNKNOWN_Z,  // 9_48 unknown 130
  UNKNOWN_Z,  // 9_49 unknown 130
  0.202974145508,  // 6c_1 good 130
  0.202974145508,  // 6c_2 good 130
  0.199505615234,  // 7c_1 good 130
  0.197113525391,  // 8c_1 good 130
  0.197113525391,  // 8c_2 good 130
  0.196236425781,  // 8c_3 good 130
  0.196236425781,  // 8c_4 good 130
  0.195518798828}; // 8c_5 good 130

double zvalue_135 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 135
  0.210284970093,  // 3_1 good 135
  0.208316479492,  // 4_1 good 135
  0.20702076416,  // 5_1 good 135
  0.206323071289,  // 5_2 good 135
  0.204030651855,  // 6_1 good 135
  0.203711706543,  // 6_2 good 135
  0.203332958984,  // 6_3 good 135
  0.202575463867,  // 7_1 good 135
  0.201578759766,  // 7_2 good 135
  0.201339550781,  // 7_3 good 135
  0.201020605469,  // 7_4 good 135
  0.200621923828,  // 7_5 good 135
  0.200462451172,  // 7_6 good 135
  0.200302978516,  // 7_7 good 135
  0.198668383789,  // 8_1 good 135
  0.198229833984,  // 8_2 good 135
  0.198389306641,  // 8_3 good 135
  0.198070361328,  // 8_4 good 135
  0.197751416016,  // 8_5 good 135
  0.197512207031,  // 8_6 good 135
  0.197751416016,  // 8_7 good 135
  0.197193261719,  // 8_8 good 135
  0.197512207031,  // 8_9 good 135
  0.197193261719,  // 8_10 good 135
  0.197193261719,  // 8_11 good 135
  0.196794580078,  // 8_12 good 135
  0.197193261719,  // 8_13 good 135
  0.196794580078,  // 8_14 good 135
  0.196674975586,  // 8_15 good 135
  0.197193261719,  // 8_16 good 135
  0.196794580078,  // 8_17 good 135
  0.196954052734,  // 8_18 good 135
  0.203851245117,  // 8_19 good 135
  0.202415991211,  // 8_20 good 135
  0.201419287109,  // 8_21 good 135
  UNKNOWN_Z,  // 9_1 unknown 135
  UNKNOWN_Z,  // 9_2 unknown 135
  UNKNOWN_Z,  // 9_3 unknown 135
  UNKNOWN_Z,  // 9_4 unknown 135
  UNKNOWN_Z,  // 9_5 unknown 135
  UNKNOWN_Z,  // 9_6 unknown 135
  UNKNOWN_Z,  // 9_7 unknown 135
  UNKNOWN_Z,  // 9_8 unknown 135
  UNKNOWN_Z,  // 9_9 unknown 135
  UNKNOWN_Z,  // 9_10 unknown 135
  UNKNOWN_Z,  // 9_11 unknown 135
  UNKNOWN_Z,  // 9_12 unknown 135
  UNKNOWN_Z,  // 9_13 unknown 135
  UNKNOWN_Z,  // 9_14 unknown 135
  UNKNOWN_Z,  // 9_15 unknown 135
  UNKNOWN_Z,  // 9_16 unknown 135
  UNKNOWN_Z,  // 9_17 unknown 135
  UNKNOWN_Z,  // 9_18 unknown 135
  UNKNOWN_Z,  // 9_19 unknown 135
  UNKNOWN_Z,  // 9_20 unknown 135
  UNKNOWN_Z,  // 9_21 unknown 135
  UNKNOWN_Z,  // 9_22 unknown 135
  UNKNOWN_Z,  // 9_23 unknown 135
  UNKNOWN_Z,  // 9_24 unknown 135
  UNKNOWN_Z,  // 9_25 unknown 135
  UNKNOWN_Z,  // 9_26 unknown 135
  UNKNOWN_Z,  // 9_27 unknown 135
  UNKNOWN_Z,  // 9_28 unknown 135
  UNKNOWN_Z,  // 9_29 unknown 135
  UNKNOWN_Z,  // 9_30 unknown 135
  UNKNOWN_Z,  // 9_31 unknown 135
  UNKNOWN_Z,  // 9_32 unknown 135
  UNKNOWN_Z,  // 9_33 unknown 135
  UNKNOWN_Z,  // 9_34 unknown 135
  UNKNOWN_Z,  // 9_35 unknown 135
  UNKNOWN_Z,  // 9_36 unknown 135
  UNKNOWN_Z,  // 9_37 unknown 135
  UNKNOWN_Z,  // 9_38 unknown 135
  UNKNOWN_Z,  // 9_39 unknown 135
  UNKNOWN_Z,  // 9_40 unknown 135
  UNKNOWN_Z,  // 9_41 unknown 135
  UNKNOWN_Z,  // 9_42 unknown 135
  UNKNOWN_Z,  // 9_43 unknown 135
  UNKNOWN_Z,  // 9_44 unknown 135
  UNKNOWN_Z,  // 9_45 unknown 135
  UNKNOWN_Z,  // 9_46 unknown 135
  UNKNOWN_Z,  // 9_47 unknown 135
  UNKNOWN_Z,  // 9_48 unknown 135
  UNKNOWN_Z,  // 9_49 unknown 135
  0.203691772461,  // 6c_1 good 135
  0.203691772461,  // 6c_2 good 135
  0.200502319336,  // 7c_1 good 135
  0.198389306641,  // 8c_1 good 135
  0.198229833984,  // 8c_2 good 135
  0.197392602539,  // 8c_3 good 135
  0.197392602539,  // 8c_4 good 135
  0.196794580078}; // 8c_5 good 135

double zvalue_140 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 140
  0.210569030762,  // 3_1 good 140
  0.208670309448,  // 4_1 good 140
  0.207359643555,  // 5_1 good 140
  0.206801489258,  // 5_2 good 140
  0.20466854248,  // 6_1 good 140
  0.204329663086,  // 6_2 good 140
  0.204130322266,  // 6_3 good 140
  0.203332958984,  // 7_1 good 140
  0.202296386719,  // 7_2 good 140
  0.202216650391,  // 7_3 good 140
  0.201937573242,  // 7_4 good 140
  0.201578759766,  // 7_5 good 140
  0.201379418945,  // 7_6 good 140
  0.201219946289,  // 7_7 good 140
  0.199744824219,  // 8_1 good 140
  0.199226538086,  // 8_2 good 140
  0.199346142578,  // 8_3 good 140
  0.199186669922,  // 8_4 good 140
  0.198787988281,  // 8_5 good 140
  0.198469042969,  // 8_6 good 140
  0.198668383789,  // 8_7 good 140
  0.198389306641,  // 8_8 good 140
  0.198668383789,  // 8_9 good 140
  0.198229833984,  // 8_10 good 140
  0.198469042969,  // 8_11 good 140
  0.198070361328,  // 8_12 good 140
  0.198229833984,  // 8_13 good 140
  0.197910888672,  // 8_14 good 140
  0.197910888672,  // 8_15 good 140
  0.198229833984,  // 8_16 good 140
  0.198030493164,  // 8_17 good 140
  0.197950756836,  // 8_18 good 140
  0.204489135742,  // 8_19 good 140
  0.203133618164,  // 8_20 good 140
  0.202216650391,  // 8_21 good 140
  UNKNOWN_Z,  // 9_1 unknown 140
  UNKNOWN_Z,  // 9_2 unknown 140
  UNKNOWN_Z,  // 9_3 unknown 140
  UNKNOWN_Z,  // 9_4 unknown 140
  UNKNOWN_Z,  // 9_5 unknown 140
  UNKNOWN_Z,  // 9_6 unknown 140
  UNKNOWN_Z,  // 9_7 unknown 140
  UNKNOWN_Z,  // 9_8 unknown 140
  UNKNOWN_Z,  // 9_9 unknown 140
  UNKNOWN_Z,  // 9_10 unknown 140
  UNKNOWN_Z,  // 9_11 unknown 140
  UNKNOWN_Z,  // 9_12 unknown 140
  UNKNOWN_Z,  // 9_13 unknown 140
  UNKNOWN_Z,  // 9_14 unknown 140
  UNKNOWN_Z,  // 9_15 unknown 140
  UNKNOWN_Z,  // 9_16 unknown 140
  UNKNOWN_Z,  // 9_17 unknown 140
  UNKNOWN_Z,  // 9_18 unknown 140
  UNKNOWN_Z,  // 9_19 unknown 140
  UNKNOWN_Z,  // 9_20 unknown 140
  UNKNOWN_Z,  // 9_21 unknown 140
  UNKNOWN_Z,  // 9_22 unknown 140
  UNKNOWN_Z,  // 9_23 unknown 140
  UNKNOWN_Z,  // 9_24 unknown 140
  UNKNOWN_Z,  // 9_25 unknown 140
  UNKNOWN_Z,  // 9_26 unknown 140
  UNKNOWN_Z,  // 9_27 unknown 140
  UNKNOWN_Z,  // 9_28 unknown 140
  UNKNOWN_Z,  // 9_29 unknown 140
  UNKNOWN_Z,  // 9_30 unknown 140
  UNKNOWN_Z,  // 9_31 unknown 140
  UNKNOWN_Z,  // 9_32 unknown 140
  UNKNOWN_Z,  // 9_33 unknown 140
  UNKNOWN_Z,  // 9_34 unknown 140
  UNKNOWN_Z,  // 9_35 unknown 140
  UNKNOWN_Z,  // 9_36 unknown 140
  UNKNOWN_Z,  // 9_37 unknown 140
  UNKNOWN_Z,  // 9_38 unknown 140
  UNKNOWN_Z,  // 9_39 unknown 140
  UNKNOWN_Z,  // 9_40 unknown 140
  UNKNOWN_Z,  // 9_41 unknown 140
  UNKNOWN_Z,  // 9_42 unknown 140
  UNKNOWN_Z,  // 9_43 unknown 140
  UNKNOWN_Z,  // 9_44 unknown 140
  UNKNOWN_Z,  // 9_45 unknown 140
  UNKNOWN_Z,  // 9_46 unknown 140
  UNKNOWN_Z,  // 9_47 unknown 140
  UNKNOWN_Z,  // 9_48 unknown 140
  UNKNOWN_Z,  // 9_49 unknown 140
  0.204349597168,  // 6c_1 good 140
  0.204329663086,  // 6c_2 good 140
  0.201419287109,  // 7c_1 good 140
  0.199306274414,  // 8c_1 good 140
  0.199306274414,  // 8c_2 good 140
  0.198469042969,  // 8c_3 good 140
  0.198469042969,  // 8c_4 good 140
  0.197910888672}; // 8c_5 good 140

double zvalue_145 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 145
  0.210711065963,  // 3_1 failed 145
  0.208954370117,  // 4_1 good 145
  0.207768292236,  // 5_1 good 145
  0.207250006104,  // 5_2 good 145
  0.205246630859,  // 6_1 good 145
  0.204987487793,  // 6_2 good 145
  0.20466854248,  // 6_3 good 145
  0.203891113281,  // 7_1 good 145
  0.203053881836,  // 7_2 good 145
  0.202854541016,  // 7_3 good 145
  0.202615332031,  // 7_4 good 145
  0.202296386719,  // 7_5 good 145
  0.202057177734,  // 7_6 good 145
  0.201897705078,  // 7_7 good 145
  0.200582055664,  // 8_1 good 145
  0.200143505859,  // 8_2 good 145
  0.200302978516,  // 8_3 good 145
  0.200143505859,  // 8_4 good 145
  0.199744824219,  // 8_5 good 145
  0.199665087891,  // 8_6 good 145
  0.199744824219,  // 8_7 good 145
  0.199346142578,  // 8_8 good 145
  0.199665087891,  // 8_9 good 145
  0.199226538086,  // 8_10 good 145
  0.199346142578,  // 8_11 good 145
  0.199027197266,  // 8_12 good 145
  0.199226538086,  // 8_13 good 145
  0.199027197266,  // 8_14 good 145
  0.198867724609,  // 8_15 good 145
  0.199306274414,  // 8_16 good 145
  0.199027197266,  // 8_17 good 145
  0.199027197266,  // 8_18 good 145
  0.205127026367,  // 8_19 good 145
  0.203851245117,  // 8_20 good 145
  0.203014013672,  // 8_21 good 145
  UNKNOWN_Z,  // 9_1 unknown 145
  UNKNOWN_Z,  // 9_2 unknown 145
  UNKNOWN_Z,  // 9_3 unknown 145
  UNKNOWN_Z,  // 9_4 unknown 145
  UNKNOWN_Z,  // 9_5 unknown 145
  UNKNOWN_Z,  // 9_6 unknown 145
  UNKNOWN_Z,  // 9_7 unknown 145
  UNKNOWN_Z,  // 9_8 unknown 145
  UNKNOWN_Z,  // 9_9 unknown 145
  UNKNOWN_Z,  // 9_10 unknown 145
  UNKNOWN_Z,  // 9_11 unknown 145
  UNKNOWN_Z,  // 9_12 unknown 145
  UNKNOWN_Z,  // 9_13 unknown 145
  UNKNOWN_Z,  // 9_14 unknown 145
  UNKNOWN_Z,  // 9_15 unknown 145
  UNKNOWN_Z,  // 9_16 unknown 145
  UNKNOWN_Z,  // 9_17 unknown 145
  UNKNOWN_Z,  // 9_18 unknown 145
  UNKNOWN_Z,  // 9_19 unknown 145
  UNKNOWN_Z,  // 9_20 unknown 145
  UNKNOWN_Z,  // 9_21 unknown 145
  UNKNOWN_Z,  // 9_22 unknown 145
  UNKNOWN_Z,  // 9_23 unknown 145
  UNKNOWN_Z,  // 9_24 unknown 145
  UNKNOWN_Z,  // 9_25 unknown 145
  UNKNOWN_Z,  // 9_26 unknown 145
  UNKNOWN_Z,  // 9_27 unknown 145
  UNKNOWN_Z,  // 9_28 unknown 145
  UNKNOWN_Z,  // 9_29 unknown 145
  UNKNOWN_Z,  // 9_30 unknown 145
  UNKNOWN_Z,  // 9_31 unknown 145
  UNKNOWN_Z,  // 9_32 unknown 145
  UNKNOWN_Z,  // 9_33 unknown 145
  UNKNOWN_Z,  // 9_34 unknown 145
  UNKNOWN_Z,  // 9_35 unknown 145
  UNKNOWN_Z,  // 9_36 unknown 145
  UNKNOWN_Z,  // 9_37 unknown 145
  UNKNOWN_Z,  // 9_38 unknown 145
  UNKNOWN_Z,  // 9_39 unknown 145
  UNKNOWN_Z,  // 9_40 unknown 145
  UNKNOWN_Z,  // 9_41 unknown 145
  UNKNOWN_Z,  // 9_42 unknown 145
  UNKNOWN_Z,  // 9_43 unknown 145
  UNKNOWN_Z,  // 9_44 unknown 145
  UNKNOWN_Z,  // 9_45 unknown 145
  UNKNOWN_Z,  // 9_46 unknown 145
  UNKNOWN_Z,  // 9_47 unknown 145
  UNKNOWN_Z,  // 9_48 unknown 145
  UNKNOWN_Z,  // 9_49 unknown 145
  0.204887817383,  // 6c_1 good 145
  0.204887817383,  // 6c_2 good 145
  0.202216650391,  // 7c_1 good 145
  0.200143505859,  // 8c_1 good 145
  0.200143505859,  // 8c_2 good 145
  0.199346142578,  // 8c_3 good 145
  0.199346142578,  // 8c_4 good 145
  0.198787988281}; // 8c_5 good 145

double zvalue_150 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 150
  0.210797265625,  // 3_1 good 150
  0.209209442139,  // 4_1 good 150
  0.208093261719,  // 5_1 good 150
  0.207609863281,  // 5_2 good 150
  0.205783691406,  // 6_1 good 150
  0.205541992187,  // 6_2 good 150
  0.205192871094,  // 6_3 good 150
  0.204575195313,  // 7_1 good 150
  0.203662109375,  // 7_2 good 150
  0.203571472168,  // 7_3 good 150
  0.203232421875,  // 7_4 good 150
  0.203044433594,  // 7_5 good 150
  0.202749023438,  // 7_6 good 150
  0.202708740234,  // 7_7 good 150
  0.201352539063,  // 8_1 good 150
  0.201083984375,  // 8_2 good 150
  0.201137695313,  // 8_3 good 150
  0.201110839844,  // 8_4 good 150
  0.200600585938,  // 8_5 good 150
  0.200520019531,  // 8_6 good 150
  0.200614013672,  // 8_7 good 150
  0.200278320313,  // 8_8 good 150
  0.200466308594,  // 8_9 good 150
  0.200224609375,  // 8_10 good 150
  0.200278320313,  // 8_11 good 150
  0.200009765625,  // 8_12 good 150
  0.200197753906,  // 8_13 good 150
  0.200009765625,  // 8_14 good 150
  0.199741210937,  // 8_15 good 150
  0.200090332031,  // 8_16 good 150
  0.199848632812,  // 8_17 good 150
  0.199794921875,  // 8_18 good 150
  0.205703138113,  // 8_19 good 150
  0.204521484375,  // 8_20 good 150
  0.203581542969,  // 8_21 good 150
  0.200063769531,  // 9_1 good 150
  0.198787988281,  // 9_2 good 150
  0.198787988281,  // 9_3 good 150
  0.198469042969,  // 9_4 good 150
  0.198469042969,  // 9_5 good 150
  0.197950756836,  // 9_6 good 150
  0.197591943359,  // 9_7 good 150
  0.197512207031,  // 9_8 good 150
  0.197751416016,  // 9_9 good 150
  0.197512207031,  // 9_10 good 150
  0.197512207031,  // 9_11 good 150
  0.197512207031,  // 9_12 good 150
  0.197512207031,  // 9_13 good 150
  0.197472338867,  // 9_14 good 150
  0.197193261719,  // 9_15 good 150
  0.197193261719,  // 9_16 good 150
  0.197193261719,  // 9_17 good 150
  0.197193261719,  // 9_18 good 150
  0.197193261719,  // 9_19 good 150
  0.197193261719,  // 9_20 good 150
  0.196954052734,  // 9_21 good 150
  0.196954052734,  // 9_22 good 150
  0.196794580078,  // 9_23 good 150
  0.196794580078,  // 9_24 good 150
  0.196674975586,  // 9_25 good 150
  0.196954052734,  // 9_26 good 150
  0.196794580078,  // 9_27 good 150
  0.196635107422,  // 9_28 good 150
  0.197113525391,  // 9_29 good 150
  0.196635107422,  // 9_30 good 150
  0.196635107422,  // 9_31 good 150
  0.196794580078,  // 9_32 good 150
  0.196674975586,  // 9_33 good 150
  0.196834448242,  // 9_34 good 150
  0.197910888672,  // 9_35 good 150
  0.197193261719,  // 9_36 good 150
  0.196954052734,  // 9_37 good 150
  0.196595239258,  // 9_38 good 150
  0.196954052734,  // 9_39 good 150
  0.197193261719,  // 9_40 good 150
  0.197193261719,  // 9_41 good 150
  0.202695068359,  // 9_42 good 150
  0.201778100586,  // 9_43 good 150
  0.201419287109,  // 9_44 good 150
  0.200781396484,  // 9_45 good 150
  0.202854541016,  // 9_46 good 150
  0.200781396484,  // 9_47 good 150
  0.200621923828,  // 9_48 good 150
  0.201020605469,  // 9_49 good 150
  0.205369125366,  // 6c_1 good 150
  0.20541151123,  // 6c_2 good 150
  0.202788574219,  // 7c_1 good 150
  0.201033300781,  // 8c_1 good 150
  0.201033300781,  // 8c_2 good 150
  0.200335180664,  // 8c_3 good 150
  0.200235449219,  // 8c_4 good 150
  0.199876416016}; // 8c_5 good 150

double zvalue_155 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 155
  0.21107735472,  // 3_1 failed 155
  0.209497573853,  // 4_1 good 155
  0.208436161852,  // 5_1 good 155
  0.207927764893,  // 5_2 good 155
  0.206243334961,  // 6_1 good 155
  0.20594432373,  // 6_2 good 155
  0.205685180664,  // 6_3 good 155
  0.205047290039,  // 7_1 good 155
  0.204329663086,  // 7_2 good 155
  0.204130322266,  // 7_3 good 155
  0.203970849609,  // 7_4 good 155
  0.203691772461,  // 7_5 good 155
  0.203492431641,  // 7_6 good 155
  0.203332958984,  // 7_7 good 155
  0.20211697998,  // 8_1 good 155
  0.201857836914,  // 8_2 good 155
  0.201937573242,  // 8_3 good 155
  0.201778100586,  // 8_4 good 155
  0.201419287109,  // 8_5 good 155
  0.201299682617,  // 8_6 good 155
  0.201339550781,  // 8_7 good 155
  0.201140209961,  // 8_8 good 155
  0.201219946289,  // 8_9 good 155
  0.201020605469,  // 8_10 good 155
  0.201060473633,  // 8_11 good 155
  0.200781396484,  // 8_12 good 155
  0.201020605469,  // 8_13 good 155
  0.200781396484,  // 8_14 good 155
  0.200621923828,  // 8_15 good 155
  0.200940869141,  // 8_16 good 155
  0.200781396484,  // 8_17 good 155
  0.200781396484,  // 8_18 good 155
  0.206043994141,  // 8_19 good 155
  0.205027355957,  // 8_20 good 155
  0.204289794922,  // 8_21 good 155
  UNKNOWN_Z,  // 9_1 unknown 155
  UNKNOWN_Z,  // 9_2 unknown 155
  UNKNOWN_Z,  // 9_3 unknown 155
  UNKNOWN_Z,  // 9_4 unknown 155
  UNKNOWN_Z,  // 9_5 unknown 155
  UNKNOWN_Z,  // 9_6 unknown 155
  UNKNOWN_Z,  // 9_7 unknown 155
  UNKNOWN_Z,  // 9_8 unknown 155
  UNKNOWN_Z,  // 9_9 unknown 155
  UNKNOWN_Z,  // 9_10 unknown 155
  UNKNOWN_Z,  // 9_11 unknown 155
  UNKNOWN_Z,  // 9_12 unknown 155
  UNKNOWN_Z,  // 9_13 unknown 155
  UNKNOWN_Z,  // 9_14 unknown 155
  UNKNOWN_Z,  // 9_15 unknown 155
  UNKNOWN_Z,  // 9_16 unknown 155
  UNKNOWN_Z,  // 9_17 unknown 155
  UNKNOWN_Z,  // 9_18 unknown 155
  UNKNOWN_Z,  // 9_19 unknown 155
  UNKNOWN_Z,  // 9_20 unknown 155
  UNKNOWN_Z,  // 9_21 unknown 155
  UNKNOWN_Z,  // 9_22 unknown 155
  UNKNOWN_Z,  // 9_23 unknown 155
  UNKNOWN_Z,  // 9_24 unknown 155
  UNKNOWN_Z,  // 9_25 unknown 155
  UNKNOWN_Z,  // 9_26 unknown 155
  UNKNOWN_Z,  // 9_27 unknown 155
  UNKNOWN_Z,  // 9_28 unknown 155
  UNKNOWN_Z,  // 9_29 unknown 155
  UNKNOWN_Z,  // 9_30 unknown 155
  UNKNOWN_Z,  // 9_31 unknown 155
  UNKNOWN_Z,  // 9_32 unknown 155
  UNKNOWN_Z,  // 9_33 unknown 155
  UNKNOWN_Z,  // 9_34 unknown 155
  UNKNOWN_Z,  // 9_35 unknown 155
  UNKNOWN_Z,  // 9_36 unknown 155
  UNKNOWN_Z,  // 9_37 unknown 155
  UNKNOWN_Z,  // 9_38 unknown 155
  UNKNOWN_Z,  // 9_39 unknown 155
  UNKNOWN_Z,  // 9_40 unknown 155
  UNKNOWN_Z,  // 9_41 unknown 155
  UNKNOWN_Z,  // 9_42 unknown 155
  UNKNOWN_Z,  // 9_43 unknown 155
  UNKNOWN_Z,  // 9_44 unknown 155
  UNKNOWN_Z,  // 9_45 unknown 155
  UNKNOWN_Z,  // 9_46 unknown 155
  UNKNOWN_Z,  // 9_47 unknown 155
  UNKNOWN_Z,  // 9_48 unknown 155
  UNKNOWN_Z,  // 9_49 unknown 155
  0.20584465332,  // 6c_1 good 155
  0.205794818115,  // 6c_2 good 155
  0.203492431641,  // 7c_1 good 155
  0.201778100586,  // 8c_1 good 155
  0.201837902832,  // 8c_2 good 155
  0.201060473633,  // 8c_3 good 155
  0.201140209961,  // 8c_4 good 155
  0.200621923828}; // 8c_5 good 155

double zvalue_160 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 160
  0.211087316895,  // 3_1 good 160
  0.209721832275,  // 4_1 good 160
  0.208635424805,  // 5_1 good 160
  0.20819438324,  // 5_2 good 160
  0.20662208252,  // 6_1 good 160
  0.206402807617,  // 6_2 good 160
  0.206083862305,  // 6_3 good 160
  0.205525708008,  // 7_1 good 160
  0.204887817383,  // 7_2 good 160
  0.204683493042,  // 7_3 good 160
  0.204489135742,  // 7_4 good 160
  0.20417019043,  // 7_5 good 160
  0.204030651855,  // 7_6 good 160
  0.203811376953,  // 7_7 good 160
  0.202854541016,  // 8_1 good 160
  0.202575463867,  // 8_2 good 160
  0.202615332031,  // 8_3 good 160
  0.202415991211,  // 8_4 good 160
  0.202216650391,  // 8_5 good 160
  0.202057177734,  // 8_6 good 160
  0.202057177734,  // 8_7 good 160
  0.201798034668,  // 8_8 good 160
  0.201897705078,  // 8_9 good 160
  0.201778100586,  // 8_10 good 160
  0.201778100586,  // 8_11 good 160
  0.201578759766,  // 8_12 good 160
  0.201698364258,  // 8_13 good 160
  0.201578759766,  // 8_14 good 160
  0.201259814453,  // 8_15 good 160
  0.201578759766,  // 8_16 good 160
  0.201479089355,  // 8_17 good 160
  0.201419287109,  // 8_18 good 160
  0.20643270874,  // 8_19 good 160
  0.205525708008,  // 8_20 good 160
  0.204768212891,  // 8_21 good 160
  UNKNOWN_Z,  // 9_1 unknown 160
  UNKNOWN_Z,  // 9_2 unknown 160
  UNKNOWN_Z,  // 9_3 unknown 160
  UNKNOWN_Z,  // 9_4 unknown 160
  UNKNOWN_Z,  // 9_5 unknown 160
  UNKNOWN_Z,  // 9_6 unknown 160
  UNKNOWN_Z,  // 9_7 unknown 160
  UNKNOWN_Z,  // 9_8 unknown 160
  UNKNOWN_Z,  // 9_9 unknown 160
  UNKNOWN_Z,  // 9_10 unknown 160
  UNKNOWN_Z,  // 9_11 unknown 160
  UNKNOWN_Z,  // 9_12 unknown 160
  UNKNOWN_Z,  // 9_13 unknown 160
  UNKNOWN_Z,  // 9_14 unknown 160
  UNKNOWN_Z,  // 9_15 unknown 160
  UNKNOWN_Z,  // 9_16 unknown 160
  UNKNOWN_Z,  // 9_17 unknown 160
  UNKNOWN_Z,  // 9_18 unknown 160
  UNKNOWN_Z,  // 9_19 unknown 160
  UNKNOWN_Z,  // 9_20 unknown 160
  UNKNOWN_Z,  // 9_21 unknown 160
  UNKNOWN_Z,  // 9_22 unknown 160
  UNKNOWN_Z,  // 9_23 unknown 160
  UNKNOWN_Z,  // 9_24 unknown 160
  UNKNOWN_Z,  // 9_25 unknown 160
  UNKNOWN_Z,  // 9_26 unknown 160
  UNKNOWN_Z,  // 9_27 unknown 160
  UNKNOWN_Z,  // 9_28 unknown 160
  UNKNOWN_Z,  // 9_29 unknown 160
  UNKNOWN_Z,  // 9_30 unknown 160
  UNKNOWN_Z,  // 9_31 unknown 160
  UNKNOWN_Z,  // 9_32 unknown 160
  UNKNOWN_Z,  // 9_33 unknown 160
  UNKNOWN_Z,  // 9_34 unknown 160
  UNKNOWN_Z,  // 9_35 unknown 160
  UNKNOWN_Z,  // 9_36 unknown 160
  UNKNOWN_Z,  // 9_37 unknown 160
  UNKNOWN_Z,  // 9_38 unknown 160
  UNKNOWN_Z,  // 9_39 unknown 160
  UNKNOWN_Z,  // 9_40 unknown 160
  UNKNOWN_Z,  // 9_41 unknown 160
  UNKNOWN_Z,  // 9_42 unknown 160
  UNKNOWN_Z,  // 9_43 unknown 160
  UNKNOWN_Z,  // 9_44 unknown 160
  UNKNOWN_Z,  // 9_45 unknown 160
  UNKNOWN_Z,  // 9_46 unknown 160
  UNKNOWN_Z,  // 9_47 unknown 160
  UNKNOWN_Z,  // 9_48 unknown 160
  UNKNOWN_Z,  // 9_49 unknown 160
  0.206273236084,  // 6c_1 good 160
  0.206263269043,  // 6c_2 good 160
  0.204130322266,  // 7c_1 good 160
  0.202415991211,  // 8c_1 good 160
  0.202415991211,  // 8c_2 good 160
  0.201778100586,  // 8c_3 good 160
  0.201778100586,  // 8c_4 good 160
  0.201419287109}; // 8c_5 good 160

double zvalue_165 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 165
  0.21126738064,  // 3_1 failed 165
  0.209911206055,  // 4_1 good 165
  0.208874633789,  // 5_1 good 165
  0.20849090271,  // 5_2 good 165
  0.206960961914,  // 6_1 good 165
  0.206801489258,  // 6_2 good 165
  0.206604640198,  // 6_3 good 165
  0.20594432373,  // 7_1 good 165
  0.205306433105,  // 7_2 good 165
  0.205206762695,  // 7_3 good 165
  0.205027355957,  // 7_4 good 165
  0.204748278809,  // 7_5 good 165
  0.204608740234,  // 7_6 good 165
  0.20456887207,  // 7_7 good 165
  0.203492431641,  // 8_1 good 165
  0.203073815918,  // 8_2 good 165
  0.203173486328,  // 8_3 good 165
  0.203053881836,  // 8_4 good 165
  0.202854541016,  // 8_5 good 165
  0.202695068359,  // 8_6 good 165
  0.202695068359,  // 8_7 good 165
  0.202495727539,  // 8_8 good 165
  0.202695068359,  // 8_9 good 165
  0.202415991211,  // 8_10 good 165
  0.202435925293,  // 8_11 good 165
  0.202256518555,  // 8_12 good 165
  0.202415991211,  // 8_13 good 165
  0.202216650391,  // 8_14 good 165
  0.202057177734,  // 8_15 good 165
  0.202296386719,  // 8_16 good 165
  0.202216650391,  // 8_17 good 165
  0.202057177734,  // 8_18 good 165
  0.20692483139,  // 8_19 good 165
  0.205884521484,  // 8_20 good 165
  0.205246630859,  // 8_21 good 165
  UNKNOWN_Z,  // 9_1 unknown 165
  UNKNOWN_Z,  // 9_2 unknown 165
  UNKNOWN_Z,  // 9_3 unknown 165
  UNKNOWN_Z,  // 9_4 unknown 165
  UNKNOWN_Z,  // 9_5 unknown 165
  UNKNOWN_Z,  // 9_6 unknown 165
  UNKNOWN_Z,  // 9_7 unknown 165
  UNKNOWN_Z,  // 9_8 unknown 165
  UNKNOWN_Z,  // 9_9 unknown 165
  UNKNOWN_Z,  // 9_10 unknown 165
  UNKNOWN_Z,  // 9_11 unknown 165
  UNKNOWN_Z,  // 9_12 unknown 165
  UNKNOWN_Z,  // 9_13 unknown 165
  UNKNOWN_Z,  // 9_14 unknown 165
  UNKNOWN_Z,  // 9_15 unknown 165
  UNKNOWN_Z,  // 9_16 unknown 165
  UNKNOWN_Z,  // 9_17 unknown 165
  UNKNOWN_Z,  // 9_18 unknown 165
  UNKNOWN_Z,  // 9_19 unknown 165
  UNKNOWN_Z,  // 9_20 unknown 165
  UNKNOWN_Z,  // 9_21 unknown 165
  UNKNOWN_Z,  // 9_22 unknown 165
  UNKNOWN_Z,  // 9_23 unknown 165
  UNKNOWN_Z,  // 9_24 unknown 165
  UNKNOWN_Z,  // 9_25 unknown 165
  UNKNOWN_Z,  // 9_26 unknown 165
  UNKNOWN_Z,  // 9_27 unknown 165
  UNKNOWN_Z,  // 9_28 unknown 165
  UNKNOWN_Z,  // 9_29 unknown 165
  UNKNOWN_Z,  // 9_30 unknown 165
  UNKNOWN_Z,  // 9_31 unknown 165
  UNKNOWN_Z,  // 9_32 unknown 165
  UNKNOWN_Z,  // 9_33 unknown 165
  UNKNOWN_Z,  // 9_34 unknown 165
  UNKNOWN_Z,  // 9_35 unknown 165
  UNKNOWN_Z,  // 9_36 unknown 165
  UNKNOWN_Z,  // 9_37 unknown 165
  UNKNOWN_Z,  // 9_38 unknown 165
  UNKNOWN_Z,  // 9_39 unknown 165
  UNKNOWN_Z,  // 9_40 unknown 165
  UNKNOWN_Z,  // 9_41 unknown 165
  UNKNOWN_Z,  // 9_42 unknown 165
  UNKNOWN_Z,  // 9_43 unknown 165
  UNKNOWN_Z,  // 9_44 unknown 165
  UNKNOWN_Z,  // 9_45 unknown 165
  UNKNOWN_Z,  // 9_46 unknown 165
  UNKNOWN_Z,  // 9_47 unknown 165
  UNKNOWN_Z,  // 9_48 unknown 165
  UNKNOWN_Z,  // 9_49 unknown 165
  0.206701818848,  // 6c_1 good 165
  0.20662208252,  // 6c_2 good 165
  0.20456887207,  // 7c_1 good 165
  0.203053881836,  // 8c_1 good 165
  0.203113684082,  // 8c_2 good 165
  0.202495727539,  // 8c_3 good 165
  0.202495727539,  // 8c_4 good 165
  0.202216650391}; // 8c_5 good 165

double zvalue_170 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 170
  0.211356427002,  // 3_1 good 170
  0.210140447998,  // 4_1 good 170
  0.20917364502,  // 5_1 good 170
  0.208715161133,  // 5_2 good 170
  0.207319775391,  // 6_1 good 170
  0.207160302734,  // 6_2 good 170
  0.206901159668,  // 6_3 good 170
  0.206402807617,  // 7_1 good 170
  0.205764916992,  // 7_2 good 170
  0.205685180664,  // 7_3 good 170
  0.205406103516,  // 7_4 good 170
  0.205087158203,  // 7_5 good 170
  0.205067224121,  // 7_6 good 170
  0.204967553711,  // 7_7 good 170
  0.203970849609,  // 8_1 good 170
  0.203691772461,  // 8_2 good 170
  0.203851245117,  // 8_3 good 170
  0.203612036133,  // 8_4 good 170
  0.203332958984,  // 8_5 good 170
  0.203173486328,  // 8_6 good 170
  0.203332958984,  // 8_7 good 170
  0.203073815918,  // 8_8 good 170
  0.203213354492,  // 8_9 good 170
  0.203113684082,  // 8_10 good 170
  0.203133618164,  // 8_11 good 170
  0.202854541016,  // 8_12 good 170
  0.203014013672,  // 8_13 good 170
  0.202854541016,  // 8_14 good 170
  0.20279473877,  // 8_15 good 170
  0.202854541016,  // 8_16 good 170
  0.202854541016,  // 8_17 good 170
  0.202695068359,  // 8_18 good 170
  0.20722010498,  // 8_19 good 170
  0.206323071289,  // 8_20 good 170
  0.205685180664,  // 8_21 good 170
  UNKNOWN_Z,  // 9_1 unknown 170
  UNKNOWN_Z,  // 9_2 unknown 170
  UNKNOWN_Z,  // 9_3 unknown 170
  UNKNOWN_Z,  // 9_4 unknown 170
  UNKNOWN_Z,  // 9_5 unknown 170
  UNKNOWN_Z,  // 9_6 unknown 170
  UNKNOWN_Z,  // 9_7 unknown 170
  UNKNOWN_Z,  // 9_8 unknown 170
  UNKNOWN_Z,  // 9_9 unknown 170
  UNKNOWN_Z,  // 9_10 unknown 170
  UNKNOWN_Z,  // 9_11 unknown 170
  UNKNOWN_Z,  // 9_12 unknown 170
  UNKNOWN_Z,  // 9_13 unknown 170
  UNKNOWN_Z,  // 9_14 unknown 170
  UNKNOWN_Z,  // 9_15 unknown 170
  UNKNOWN_Z,  // 9_16 unknown 170
  UNKNOWN_Z,  // 9_17 unknown 170
  UNKNOWN_Z,  // 9_18 unknown 170
  UNKNOWN_Z,  // 9_19 unknown 170
  UNKNOWN_Z,  // 9_20 unknown 170
  UNKNOWN_Z,  // 9_21 unknown 170
  UNKNOWN_Z,  // 9_22 unknown 170
  UNKNOWN_Z,  // 9_23 unknown 170
  UNKNOWN_Z,  // 9_24 unknown 170
  UNKNOWN_Z,  // 9_25 unknown 170
  UNKNOWN_Z,  // 9_26 unknown 170
  UNKNOWN_Z,  // 9_27 unknown 170
  UNKNOWN_Z,  // 9_28 unknown 170
  UNKNOWN_Z,  // 9_29 unknown 170
  UNKNOWN_Z,  // 9_30 unknown 170
  UNKNOWN_Z,  // 9_31 unknown 170
  UNKNOWN_Z,  // 9_32 unknown 170
  UNKNOWN_Z,  // 9_33 unknown 170
  UNKNOWN_Z,  // 9_34 unknown 170
  UNKNOWN_Z,  // 9_35 unknown 170
  UNKNOWN_Z,  // 9_36 unknown 170
  UNKNOWN_Z,  // 9_37 unknown 170
  UNKNOWN_Z,  // 9_38 unknown 170
  UNKNOWN_Z,  // 9_39 unknown 170
  UNKNOWN_Z,  // 9_40 unknown 170
  UNKNOWN_Z,  // 9_41 unknown 170
  UNKNOWN_Z,  // 9_42 unknown 170
  UNKNOWN_Z,  // 9_43 unknown 170
  UNKNOWN_Z,  // 9_44 unknown 170
  UNKNOWN_Z,  // 9_45 unknown 170
  UNKNOWN_Z,  // 9_46 unknown 170
  UNKNOWN_Z,  // 9_47 unknown 170
  UNKNOWN_Z,  // 9_48 unknown 170
  UNKNOWN_Z,  // 9_49 unknown 170
  0.207040698242,  // 6c_1 good 170
  0.206960961914,  // 6c_2 good 170
  0.205027355957,  // 7c_1 good 170
  0.203612036133,  // 8c_1 good 170
  0.203572167969,  // 8c_2 good 170
  0.203053881836,  // 8c_3 good 170
  0.203113684082,  // 8c_4 good 170
  0.202695068359}; // 8c_5 good 170

double zvalue_175 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 175
  0.211404704857,  // 3_1 good 175
  0.210301111317,  // 4_1 good 175
  0.209415011597,  // 5_1 good 175
  0.209039526367,  // 5_2 good 175
  0.207613720703,  // 6_1 good 175
  0.207475292969,  // 6_2 good 175
  0.207281494141,  // 6_3 good 175
  0.206672412109,  // 7_1 good 175
  0.206146386719,  // 7_2 good 175
  0.205980273438,  // 7_3 good 175
  0.205841845703,  // 7_4 good 175
  0.205592675781,  // 7_5 good 175
  0.205481933594,  // 7_6 good 175
  0.205398876953,  // 7_7 good 175
  0.204485253906,  // 8_1 good 175
  0.204208398437,  // 8_2 good 175
  0.204374511719,  // 8_3 good 175
  0.204125341797,  // 8_4 good 175
  0.203931542969,  // 8_5 good 175
  0.203820800781,  // 8_6 good 175
  0.203931542969,  // 8_7 good 175
  0.203682373047,  // 8_8 good 175
  0.203820800781,  // 8_9 good 175
  0.203599316406,  // 8_10 good 175
  0.203627001953,  // 8_11 good 175
  0.203488574219,  // 8_12 good 175
  0.203488574219,  // 8_13 good 175
  0.203322460938,  // 8_14 good 175
  0.203322460938,  // 8_15 good 175
  0.203488574219,  // 8_16 good 175
  0.203322460938,  // 8_17 good 175
  0.203322460938,  // 8_18 good 175
  0.207502978516,  // 8_19 good 175
  0.206672412109,  // 8_20 good 175
  0.206104858398,  // 8_21 good 175
  UNKNOWN_Z,  // 9_1 unknown 175
  UNKNOWN_Z,  // 9_2 unknown 175
  UNKNOWN_Z,  // 9_3 unknown 175
  UNKNOWN_Z,  // 9_4 unknown 175
  UNKNOWN_Z,  // 9_5 unknown 175
  UNKNOWN_Z,  // 9_6 unknown 175
  UNKNOWN_Z,  // 9_7 unknown 175
  UNKNOWN_Z,  // 9_8 unknown 175
  UNKNOWN_Z,  // 9_9 unknown 175
  UNKNOWN_Z,  // 9_10 unknown 175
  UNKNOWN_Z,  // 9_11 unknown 175
  UNKNOWN_Z,  // 9_12 unknown 175
  UNKNOWN_Z,  // 9_13 unknown 175
  UNKNOWN_Z,  // 9_14 unknown 175
  UNKNOWN_Z,  // 9_15 unknown 175
  UNKNOWN_Z,  // 9_16 unknown 175
  UNKNOWN_Z,  // 9_17 unknown 175
  UNKNOWN_Z,  // 9_18 unknown 175
  UNKNOWN_Z,  // 9_19 unknown 175
  UNKNOWN_Z,  // 9_20 unknown 175
  UNKNOWN_Z,  // 9_21 unknown 175
  UNKNOWN_Z,  // 9_22 unknown 175
  UNKNOWN_Z,  // 9_23 unknown 175
  UNKNOWN_Z,  // 9_24 unknown 175
  UNKNOWN_Z,  // 9_25 unknown 175
  UNKNOWN_Z,  // 9_26 unknown 175
  UNKNOWN_Z,  // 9_27 unknown 175
  UNKNOWN_Z,  // 9_28 unknown 175
  UNKNOWN_Z,  // 9_29 unknown 175
  UNKNOWN_Z,  // 9_30 unknown 175
  UNKNOWN_Z,  // 9_31 unknown 175
  UNKNOWN_Z,  // 9_32 unknown 175
  UNKNOWN_Z,  // 9_33 unknown 175
  UNKNOWN_Z,  // 9_34 unknown 175
  UNKNOWN_Z,  // 9_35 unknown 175
  UNKNOWN_Z,  // 9_36 unknown 175
  UNKNOWN_Z,  // 9_37 unknown 175
  UNKNOWN_Z,  // 9_38 unknown 175
  UNKNOWN_Z,  // 9_39 unknown 175
  UNKNOWN_Z,  // 9_40 unknown 175
  UNKNOWN_Z,  // 9_41 unknown 175
  UNKNOWN_Z,  // 9_42 unknown 175
  UNKNOWN_Z,  // 9_43 unknown 175
  UNKNOWN_Z,  // 9_44 unknown 175
  UNKNOWN_Z,  // 9_45 unknown 175
  UNKNOWN_Z,  // 9_46 unknown 175
  UNKNOWN_Z,  // 9_47 unknown 175
  UNKNOWN_Z,  // 9_48 unknown 175
  UNKNOWN_Z,  // 9_49 unknown 175
  0.207339709473,  // 6c_1 good 175
  0.207319775391,  // 6c_2 good 175
  0.205406103516,  // 7c_1 good 175
  0.204130322266,  // 8c_1 good 175
  0.204130322266,  // 8c_2 good 175
  0.203651904297,  // 8c_3 good 175
  0.203691772461,  // 8c_4 good 175
  0.203213354492}; // 8c_5 good 175

double zvalue_180 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 180
  0.211630603364,  // 3_1 failed 180
  0.210459393311,  // 4_1 good 180
  0.20959358449,  // 5_1 good 180
  0.209238430786,  // 5_2 good 180
  0.207937731934,  // 6_1 good 180
  0.207758325195,  // 6_2 good 180
  0.207578918457,  // 6_3 good 180
  0.207110467529,  // 7_1 good 180
  0.206562280273,  // 7_2 good 180
  0.206402807617,  // 7_3 good 180
  0.206183532715,  // 7_4 good 180
  0.206043994141,  // 7_5 good 180
  0.205864587402,  // 7_6 good 180
  0.205764916992,  // 7_7 good 180
  0.205027355957,  // 8_1 good 180
  0.204708410645,  // 8_2 good 180
  0.204768212891,  // 8_3 good 180
  0.204608740234,  // 8_4 good 180
  0.204389465332,  // 8_5 good 180
  0.204249926758,  // 8_6 good 180
  0.204329663086,  // 8_7 good 180
  0.204130322266,  // 8_8 good 180
  0.204249926758,  // 8_9 good 180
  0.204130322266,  // 8_10 good 180
  0.204130322266,  // 8_11 good 180
  0.203970849609,  // 8_12 good 180
  0.204030651855,  // 8_13 good 180
  0.203970849609,  // 8_14 good 180
  0.203851245117,  // 8_15 good 180
  0.203970849609,  // 8_16 good 180
  0.203930981445,  // 8_17 good 180
  0.203851245117,  // 8_18 good 180
  0.20789786377,  // 8_19 good 180
  0.207085549927,  // 8_20 good 180
  0.206482543945,  // 8_21 good 180
  UNKNOWN_Z,  // 9_1 unknown 180
  UNKNOWN_Z,  // 9_2 unknown 180
  UNKNOWN_Z,  // 9_3 unknown 180
  UNKNOWN_Z,  // 9_4 unknown 180
  UNKNOWN_Z,  // 9_5 unknown 180
  UNKNOWN_Z,  // 9_6 unknown 180
  UNKNOWN_Z,  // 9_7 unknown 180
  UNKNOWN_Z,  // 9_8 unknown 180
  UNKNOWN_Z,  // 9_9 unknown 180
  UNKNOWN_Z,  // 9_10 unknown 180
  UNKNOWN_Z,  // 9_11 unknown 180
  UNKNOWN_Z,  // 9_12 unknown 180
  UNKNOWN_Z,  // 9_13 unknown 180
  UNKNOWN_Z,  // 9_14 unknown 180
  UNKNOWN_Z,  // 9_15 unknown 180
  UNKNOWN_Z,  // 9_16 unknown 180
  UNKNOWN_Z,  // 9_17 unknown 180
  UNKNOWN_Z,  // 9_18 unknown 180
  UNKNOWN_Z,  // 9_19 unknown 180
  UNKNOWN_Z,  // 9_20 unknown 180
  UNKNOWN_Z,  // 9_21 unknown 180
  UNKNOWN_Z,  // 9_22 unknown 180
  UNKNOWN_Z,  // 9_23 unknown 180
  UNKNOWN_Z,  // 9_24 unknown 180
  UNKNOWN_Z,  // 9_25 unknown 180
  UNKNOWN_Z,  // 9_26 unknown 180
  UNKNOWN_Z,  // 9_27 unknown 180
  UNKNOWN_Z,  // 9_28 unknown 180
  UNKNOWN_Z,  // 9_29 unknown 180
  UNKNOWN_Z,  // 9_30 unknown 180
  UNKNOWN_Z,  // 9_31 unknown 180
  UNKNOWN_Z,  // 9_32 unknown 180
  UNKNOWN_Z,  // 9_33 unknown 180
  UNKNOWN_Z,  // 9_34 unknown 180
  UNKNOWN_Z,  // 9_35 unknown 180
  UNKNOWN_Z,  // 9_36 unknown 180
  UNKNOWN_Z,  // 9_37 unknown 180
  UNKNOWN_Z,  // 9_38 unknown 180
  UNKNOWN_Z,  // 9_39 unknown 180
  UNKNOWN_Z,  // 9_40 unknown 180
  UNKNOWN_Z,  // 9_41 unknown 180
  UNKNOWN_Z,  // 9_42 unknown 180
  UNKNOWN_Z,  // 9_43 unknown 180
  UNKNOWN_Z,  // 9_44 unknown 180
  UNKNOWN_Z,  // 9_45 unknown 180
  UNKNOWN_Z,  // 9_46 unknown 180
  UNKNOWN_Z,  // 9_47 unknown 180
  UNKNOWN_Z,  // 9_48 unknown 180
  UNKNOWN_Z,  // 9_49 unknown 180
  0.207578918457,  // 6c_1 good 180
  0.207552754974,  // 6c_2 good 180
  0.205884521484,  // 7c_1 good 180
  0.20456887207,  // 8c_1 good 180
  0.204608740234,  // 8c_2 good 180
  0.204130322266,  // 8c_3 good 180
  0.204030651855,  // 8c_4 good 180
  0.203771508789}; // 8c_5 good 180

double zvalue_185 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 185
  0.211660421753,  // 3_1 good 185
  0.210591456604,  // 4_1 good 185
  0.209741766357,  // 5_1 good 185
  0.209462689209,  // 5_2 good 185
  0.208216809082,  // 6_1 good 185
  0.207977600098,  // 6_2 good 185
  0.207857995605,  // 6_3 good 185
  0.207279907227,  // 7_1 good 185
  0.206881225586,  // 7_2 good 185
  0.206746670532,  // 7_3 good 185
  0.206522412109,  // 7_4 good 185
  0.206452642822,  // 7_5 good 185
  0.206243334961,  // 7_6 good 185
  0.206163598633,  // 7_7 good 185
  0.205406103516,  // 8_1 good 185
  0.205206762695,  // 8_2 good 185
  0.205127026367,  // 8_3 good 185
  0.204987487793,  // 8_4 good 185
  0.204768212891,  // 8_5 good 185
  0.204768212891,  // 8_6 good 185
  0.204768212891,  // 8_7 good 185
  0.204608740234,  // 8_8 good 185
  0.204768212891,  // 8_9 good 185
  0.204608740234,  // 8_10 good 185
  0.204608740234,  // 8_11 good 185
  0.204409399414,  // 8_12 good 185
  0.204608740234,  // 8_13 good 185
  0.204329663086,  // 8_14 good 185
  0.204329663086,  // 8_15 good 185
  0.204489135742,  // 8_16 good 185
  0.204409399414,  // 8_17 good 185
  0.204249926758,  // 8_18 good 185
  0.20809720459,  // 8_19 good 185
  0.20731479187,  // 8_20 good 185
  0.206881225586,  // 8_21 good 185
  UNKNOWN_Z,  // 9_1 unknown 185
  UNKNOWN_Z,  // 9_2 unknown 185
  UNKNOWN_Z,  // 9_3 unknown 185
  UNKNOWN_Z,  // 9_4 unknown 185
  UNKNOWN_Z,  // 9_5 unknown 185
  UNKNOWN_Z,  // 9_6 unknown 185
  UNKNOWN_Z,  // 9_7 unknown 185
  UNKNOWN_Z,  // 9_8 unknown 185
  UNKNOWN_Z,  // 9_9 unknown 185
  UNKNOWN_Z,  // 9_10 unknown 185
  UNKNOWN_Z,  // 9_11 unknown 185
  UNKNOWN_Z,  // 9_12 unknown 185
  UNKNOWN_Z,  // 9_13 unknown 185
  UNKNOWN_Z,  // 9_14 unknown 185
  UNKNOWN_Z,  // 9_15 unknown 185
  UNKNOWN_Z,  // 9_16 unknown 185
  UNKNOWN_Z,  // 9_17 unknown 185
  UNKNOWN_Z,  // 9_18 unknown 185
  UNKNOWN_Z,  // 9_19 unknown 185
  UNKNOWN_Z,  // 9_20 unknown 185
  UNKNOWN_Z,  // 9_21 unknown 185
  UNKNOWN_Z,  // 9_22 unknown 185
  UNKNOWN_Z,  // 9_23 unknown 185
  UNKNOWN_Z,  // 9_24 unknown 185
  UNKNOWN_Z,  // 9_25 unknown 185
  UNKNOWN_Z,  // 9_26 unknown 185
  UNKNOWN_Z,  // 9_27 unknown 185
  UNKNOWN_Z,  // 9_28 unknown 185
  UNKNOWN_Z,  // 9_29 unknown 185
  UNKNOWN_Z,  // 9_30 unknown 185
  UNKNOWN_Z,  // 9_31 unknown 185
  UNKNOWN_Z,  // 9_32 unknown 185
  UNKNOWN_Z,  // 9_33 unknown 185
  UNKNOWN_Z,  // 9_34 unknown 185
  UNKNOWN_Z,  // 9_35 unknown 185
  UNKNOWN_Z,  // 9_36 unknown 185
  UNKNOWN_Z,  // 9_37 unknown 185
  UNKNOWN_Z,  // 9_38 unknown 185
  UNKNOWN_Z,  // 9_39 unknown 185
  UNKNOWN_Z,  // 9_40 unknown 185
  UNKNOWN_Z,  // 9_41 unknown 185
  UNKNOWN_Z,  // 9_42 unknown 185
  UNKNOWN_Z,  // 9_43 unknown 185
  UNKNOWN_Z,  // 9_44 unknown 185
  UNKNOWN_Z,  // 9_45 unknown 185
  UNKNOWN_Z,  // 9_46 unknown 185
  UNKNOWN_Z,  // 9_47 unknown 185
  UNKNOWN_Z,  // 9_48 unknown 185
  UNKNOWN_Z,  // 9_49 unknown 185
  0.207887896729,  // 6c_1 good 185
  0.207829340363,  // 6c_2 good 185
  0.206163598633,  // 7c_1 good 185
  0.204922702026,  // 8c_1 good 185
  0.204987487793,  // 8c_2 good 185
  0.20456887207,  // 8c_3 good 185
  0.204608740234,  // 8c_4 good 185
  0.204249926758}; // 8c_5 good 185

double zvalue_190 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 190
  0.211754490569,  // 3_1 failed 190
  0.21065375061,  // 4_1 good 190
  0.209908558559,  // 5_1 good 190
  0.209657046509,  // 5_2 good 190
  0.208436083984,  // 6_1 good 190
  0.208316479492,  // 6_2 good 190
  0.208077270508,  // 6_3 good 190
  0.207678588867,  // 7_1 good 190
  0.207160302734,  // 7_2 good 190
  0.207040698242,  // 7_3 good 190
  0.206941027832,  // 7_4 good 190
  0.206711785889,  // 7_5 good 190
  0.206612115479,  // 7_6 good 190
  0.206502478027,  // 7_7 good 190
  0.205764916992,  // 8_1 good 190
  0.205525708008,  // 8_2 good 190
  0.205605444336,  // 8_3 good 190
  0.205406103516,  // 8_4 good 190
  0.205206762695,  // 8_5 good 190
  0.205127026367,  // 8_6 good 190
  0.205246630859,  // 8_7 good 190
  0.205047290039,  // 8_8 good 190
  0.205206762695,  // 8_9 good 190
  0.204967553711,  // 8_10 good 190
  0.205047290039,  // 8_11 good 190
  0.204887817383,  // 8_12 good 190
  0.204967553711,  // 8_13 good 190
  0.204967553711,  // 8_14 good 190
  0.204768212891,  // 8_15 good 190
  0.204887817383,  // 8_16 good 190
  0.204768212891,  // 8_17 good 190
  0.204748278809,  // 8_18 good 190
  0.208316479492,  // 8_19 good 190
  0.207678588867,  // 8_20 good 190
  0.207160302734,  // 8_21 good 190
  UNKNOWN_Z,  // 9_1 unknown 190
  UNKNOWN_Z,  // 9_2 unknown 190
  UNKNOWN_Z,  // 9_3 unknown 190
  UNKNOWN_Z,  // 9_4 unknown 190
  UNKNOWN_Z,  // 9_5 unknown 190
  UNKNOWN_Z,  // 9_6 unknown 190
  UNKNOWN_Z,  // 9_7 unknown 190
  UNKNOWN_Z,  // 9_8 unknown 190
  UNKNOWN_Z,  // 9_9 unknown 190
  UNKNOWN_Z,  // 9_10 unknown 190
  UNKNOWN_Z,  // 9_11 unknown 190
  UNKNOWN_Z,  // 9_12 unknown 190
  UNKNOWN_Z,  // 9_13 unknown 190
  UNKNOWN_Z,  // 9_14 unknown 190
  UNKNOWN_Z,  // 9_15 unknown 190
  UNKNOWN_Z,  // 9_16 unknown 190
  UNKNOWN_Z,  // 9_17 unknown 190
  UNKNOWN_Z,  // 9_18 unknown 190
  UNKNOWN_Z,  // 9_19 unknown 190
  UNKNOWN_Z,  // 9_20 unknown 190
  UNKNOWN_Z,  // 9_21 unknown 190
  UNKNOWN_Z,  // 9_22 unknown 190
  UNKNOWN_Z,  // 9_23 unknown 190
  UNKNOWN_Z,  // 9_24 unknown 190
  UNKNOWN_Z,  // 9_25 unknown 190
  UNKNOWN_Z,  // 9_26 unknown 190
  UNKNOWN_Z,  // 9_27 unknown 190
  UNKNOWN_Z,  // 9_28 unknown 190
  UNKNOWN_Z,  // 9_29 unknown 190
  UNKNOWN_Z,  // 9_30 unknown 190
  UNKNOWN_Z,  // 9_31 unknown 190
  UNKNOWN_Z,  // 9_32 unknown 190
  UNKNOWN_Z,  // 9_33 unknown 190
  UNKNOWN_Z,  // 9_34 unknown 190
  UNKNOWN_Z,  // 9_35 unknown 190
  UNKNOWN_Z,  // 9_36 unknown 190
  UNKNOWN_Z,  // 9_37 unknown 190
  UNKNOWN_Z,  // 9_38 unknown 190
  UNKNOWN_Z,  // 9_39 unknown 190
  UNKNOWN_Z,  // 9_40 unknown 190
  UNKNOWN_Z,  // 9_41 unknown 190
  UNKNOWN_Z,  // 9_42 unknown 190
  UNKNOWN_Z,  // 9_43 unknown 190
  UNKNOWN_Z,  // 9_44 unknown 190
  UNKNOWN_Z,  // 9_45 unknown 190
  UNKNOWN_Z,  // 9_46 unknown 190
  UNKNOWN_Z,  // 9_47 unknown 190
  UNKNOWN_Z,  // 9_48 unknown 190
  UNKNOWN_Z,  // 9_49 unknown 190
  0.208196884733,  // 6c_1 good 190
  0.208077270508,  // 6c_2 good 190
  0.206607131958,  // 7c_1 good 190
  0.205486462784,  // 8c_1 good 190
  0.205386169434,  // 8c_2 good 190
  0.204967553711,  // 8c_3 good 190
  0.204967553711,  // 8c_4 good 190
  0.204768212891}; // 8c_5 good 190

double zvalue_195 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 195
  0.211837341598,  // 3_1 failed 195
  0.210808239746,  // 4_1 good 195
  0.210090612793,  // 5_1 good 195
  0.209811924982,  // 5_2 good 195
  0.208595556641,  // 6_1 good 195
  0.20849588623,  // 6_2 good 195
  0.208366314697,  // 6_3 good 195
  0.20789786377,  // 7_1 good 195
  0.207561476135,  // 7_2 good 195
  0.207359643555,  // 7_3 good 195
  0.20722010498,  // 7_4 good 195
  0.20712043457,  // 7_5 good 195
  0.206941027832,  // 7_6 good 195
  0.206881225586,  // 7_7 good 195
  0.206153631592,  // 8_1 good 195
  0.205884521484,  // 8_2 good 195
  0.20594432373,  // 8_3 good 195
  0.20584465332,  // 8_4 good 195
  0.205605444336,  // 8_5 good 195
  0.205525708008,  // 8_6 good 195
  0.205625378418,  // 8_7 good 195
  0.205565576172,  // 8_8 good 195
  0.205605444336,  // 8_9 good 195
  0.205386169434,  // 8_10 good 195
  0.20544597168,  // 8_11 good 195
  0.205306433105,  // 8_12 good 195
  0.205406103516,  // 8_13 good 195
  0.205246630859,  // 8_14 good 195
  0.205127026367,  // 8_15 good 195
  0.205246630859,  // 8_16 good 195
  0.205127026367,  // 8_17 good 195
  0.205127026367,  // 8_18 good 195
  0.208575622559,  // 8_19 good 195
  0.20789786377,  // 8_20 good 195
  0.207439379883,  // 8_21 good 195
  UNKNOWN_Z,  // 9_1 unknown 195
  UNKNOWN_Z,  // 9_2 unknown 195
  UNKNOWN_Z,  // 9_3 unknown 195
  UNKNOWN_Z,  // 9_4 unknown 195
  UNKNOWN_Z,  // 9_5 unknown 195
  UNKNOWN_Z,  // 9_6 unknown 195
  UNKNOWN_Z,  // 9_7 unknown 195
  UNKNOWN_Z,  // 9_8 unknown 195
  UNKNOWN_Z,  // 9_9 unknown 195
  UNKNOWN_Z,  // 9_10 unknown 195
  UNKNOWN_Z,  // 9_11 unknown 195
  UNKNOWN_Z,  // 9_12 unknown 195
  UNKNOWN_Z,  // 9_13 unknown 195
  UNKNOWN_Z,  // 9_14 unknown 195
  UNKNOWN_Z,  // 9_15 unknown 195
  UNKNOWN_Z,  // 9_16 unknown 195
  UNKNOWN_Z,  // 9_17 unknown 195
  UNKNOWN_Z,  // 9_18 unknown 195
  UNKNOWN_Z,  // 9_19 unknown 195
  UNKNOWN_Z,  // 9_20 unknown 195
  UNKNOWN_Z,  // 9_21 unknown 195
  UNKNOWN_Z,  // 9_22 unknown 195
  UNKNOWN_Z,  // 9_23 unknown 195
  UNKNOWN_Z,  // 9_24 unknown 195
  UNKNOWN_Z,  // 9_25 unknown 195
  UNKNOWN_Z,  // 9_26 unknown 195
  UNKNOWN_Z,  // 9_27 unknown 195
  UNKNOWN_Z,  // 9_28 unknown 195
  UNKNOWN_Z,  // 9_29 unknown 195
  UNKNOWN_Z,  // 9_30 unknown 195
  UNKNOWN_Z,  // 9_31 unknown 195
  UNKNOWN_Z,  // 9_32 unknown 195
  UNKNOWN_Z,  // 9_33 unknown 195
  UNKNOWN_Z,  // 9_34 unknown 195
  UNKNOWN_Z,  // 9_35 unknown 195
  UNKNOWN_Z,  // 9_36 unknown 195
  UNKNOWN_Z,  // 9_37 unknown 195
  UNKNOWN_Z,  // 9_38 unknown 195
  UNKNOWN_Z,  // 9_39 unknown 195
  UNKNOWN_Z,  // 9_40 unknown 195
  UNKNOWN_Z,  // 9_41 unknown 195
  UNKNOWN_Z,  // 9_42 unknown 195
  UNKNOWN_Z,  // 9_43 unknown 195
  UNKNOWN_Z,  // 9_44 unknown 195
  UNKNOWN_Z,  // 9_45 unknown 195
  UNKNOWN_Z,  // 9_46 unknown 195
  UNKNOWN_Z,  // 9_47 unknown 195
  UNKNOWN_Z,  // 9_48 unknown 195
  UNKNOWN_Z,  // 9_49 unknown 195
  0.208316479492,  // 6c_1 good 195
  0.208335479164,  // 6c_2 good 195
  0.206901159668,  // 7c_1 good 195
  0.205700131226,  // 8c_1 good 195
  0.205685180664,  // 8c_2 good 195
  0.205406103516,  // 8c_3 good 195
  0.20534630127,  // 8c_4 good 195
  0.205127026367}; // 8c_5 good 195

double zvalue_200 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 200
  0.211841619444,  // 3_1 good 200
  0.21094777832,  // 4_1 good 200
  0.210230151367,  // 5_1 good 200
  0.209968738556,  // 5_2 good 200
  0.208885498047,  // 6_1 good 200
  0.208724365234,  // 6_2 good 200
  0.208657226562,  // 6_3 good 200
  0.208146972656,  // 7_1 good 200
  0.207690429687,  // 7_2 good 200
  0.207704696655,  // 7_3 good 200
  0.207502441406,  // 7_4 good 200
  0.207341308594,  // 7_5 good 200
  0.207153320312,  // 7_6 good 200
  0.207153320312,  // 7_7 good 200
  0.206508789062,  // 8_1 good 200
  0.206240234375,  // 8_2 good 200
  0.206307373047,  // 8_3 good 200
  0.206213378906,  // 8_4 good 200
  0.205971679687,  // 8_5 good 200
  0.205877685547,  // 8_6 good 200
  0.205971679687,  // 8_7 good 200
  0.205756835937,  // 8_8 good 200
  0.205971679687,  // 8_9 good 200
  0.205729980469,  // 8_10 good 200
  0.205756835937,  // 8_11 good 200
  0.205676269531,  // 8_12 good 200
  0.205676269531,  // 8_13 good 200
  0.205541992187,  // 8_14 good 200
  0.205541992187,  // 8_15 good 200
  0.205649414062,  // 8_16 good 200
  0.205515136719,  // 8_17 good 200
  0.205541992187,  // 8_18 good 200
  0.208845214844,  // 8_19 good 200
  0.208200683594,  // 8_20 good 200
  0.207690429687,  // 8_21 good 200
  UNKNOWN_Z,  // 9_1 unknown 200
  UNKNOWN_Z,  // 9_2 unknown 200
  UNKNOWN_Z,  // 9_3 unknown 200
  UNKNOWN_Z,  // 9_4 unknown 200
  UNKNOWN_Z,  // 9_5 unknown 200
  UNKNOWN_Z,  // 9_6 unknown 200
  UNKNOWN_Z,  // 9_7 unknown 200
  UNKNOWN_Z,  // 9_8 unknown 200
  UNKNOWN_Z,  // 9_9 unknown 200
  UNKNOWN_Z,  // 9_10 unknown 200
  UNKNOWN_Z,  // 9_11 unknown 200
  UNKNOWN_Z,  // 9_12 unknown 200
  UNKNOWN_Z,  // 9_13 unknown 200
  UNKNOWN_Z,  // 9_14 unknown 200
  UNKNOWN_Z,  // 9_15 unknown 200
  UNKNOWN_Z,  // 9_16 unknown 200
  UNKNOWN_Z,  // 9_17 unknown 200
  UNKNOWN_Z,  // 9_18 unknown 200
  UNKNOWN_Z,  // 9_19 unknown 200
  UNKNOWN_Z,  // 9_20 unknown 200
  UNKNOWN_Z,  // 9_21 unknown 200
  UNKNOWN_Z,  // 9_22 unknown 200
  UNKNOWN_Z,  // 9_23 unknown 200
  UNKNOWN_Z,  // 9_24 unknown 200
  UNKNOWN_Z,  // 9_25 unknown 200
  UNKNOWN_Z,  // 9_26 unknown 200
  UNKNOWN_Z,  // 9_27 unknown 200
  UNKNOWN_Z,  // 9_28 unknown 200
  UNKNOWN_Z,  // 9_29 unknown 200
  UNKNOWN_Z,  // 9_30 unknown 200
  UNKNOWN_Z,  // 9_31 unknown 200
  UNKNOWN_Z,  // 9_32 unknown 200
  UNKNOWN_Z,  // 9_33 unknown 200
  UNKNOWN_Z,  // 9_34 unknown 200
  UNKNOWN_Z,  // 9_35 unknown 200
  UNKNOWN_Z,  // 9_36 unknown 200
  UNKNOWN_Z,  // 9_37 unknown 200
  UNKNOWN_Z,  // 9_38 unknown 200
  UNKNOWN_Z,  // 9_39 unknown 200
  UNKNOWN_Z,  // 9_40 unknown 200
  UNKNOWN_Z,  // 9_41 unknown 200
  UNKNOWN_Z,  // 9_42 unknown 200
  UNKNOWN_Z,  // 9_43 unknown 200
  UNKNOWN_Z,  // 9_44 unknown 200
  UNKNOWN_Z,  // 9_45 unknown 200
  UNKNOWN_Z,  // 9_46 unknown 200
  UNKNOWN_Z,  // 9_47 unknown 200
  UNKNOWN_Z,  // 9_48 unknown 200
  UNKNOWN_Z,  // 9_49 unknown 200
  0.208560671997,  // 6c_1 good 200
  0.208535754395,  // 6c_2 good 200
  0.207180236816,  // 7c_1 good 200
  0.206083862305,  // 8c_1 good 200
  0.206043994141,  // 8c_2 good 200
  0.205685180664,  // 8c_3 good 200
  0.205764916992,  // 8c_4 good 200
  0.205525708008}; // 8c_5 good 200

double zvalue_300 [NUM_KNOTS] = {
  UNKNOWN_Z,  // 0_1 unknown 300
  0.21258205671,  // 3_1 failed 300
  0.212185642961,  // 4_1 failed 300
  0.211782493672,  // 5_1 failed 300
  UNKNOWN_Z,  // 5_2 unknown 300
  UNKNOWN_Z,  // 6_1 unknown 300
  UNKNOWN_Z,  // 6_2 unknown 300
  UNKNOWN_Z,  // 6_3 unknown 300
  UNKNOWN_Z,  // 7_1 unknown 300
  UNKNOWN_Z,  // 7_2 unknown 300
  UNKNOWN_Z,  // 7_3 unknown 300
  UNKNOWN_Z,  // 7_4 unknown 300
  UNKNOWN_Z,  // 7_5 unknown 300
  UNKNOWN_Z,  // 7_6 unknown 300
  UNKNOWN_Z,  // 7_7 unknown 300
  UNKNOWN_Z,  // 8_1 unknown 300
  UNKNOWN_Z,  // 8_2 unknown 300
  UNKNOWN_Z,  // 8_3 unknown 300
  UNKNOWN_Z,  // 8_4 unknown 300
  UNKNOWN_Z,  // 8_5 unknown 300
  UNKNOWN_Z,  // 8_6 unknown 300
  UNKNOWN_Z,  // 8_7 unknown 300
  UNKNOWN_Z,  // 8_8 unknown 300
  UNKNOWN_Z,  // 8_9 unknown 300
  UNKNOWN_Z,  // 8_10 unknown 300
  UNKNOWN_Z,  // 8_11 unknown 300
  UNKNOWN_Z,  // 8_12 unknown 300
  UNKNOWN_Z,  // 8_13 unknown 300
  UNKNOWN_Z,  // 8_14 unknown 300
  UNKNOWN_Z,  // 8_15 unknown 300
  UNKNOWN_Z,  // 8_16 unknown 300
  UNKNOWN_Z,  // 8_17 unknown 300
  UNKNOWN_Z,  // 8_18 unknown 300
  UNKNOWN_Z,  // 8_19 unknown 300
  UNKNOWN_Z,  // 8_20 unknown 300
  UNKNOWN_Z,  // 8_21 unknown 300
  UNKNOWN_Z,  // 9_1 unknown 300
  UNKNOWN_Z,  // 9_2 unknown 300
  UNKNOWN_Z,  // 9_3 unknown 300
  UNKNOWN_Z,  // 9_4 unknown 300
  UNKNOWN_Z,  // 9_5 unknown 300
  UNKNOWN_Z,  // 9_6 unknown 300
  UNKNOWN_Z,  // 9_7 unknown 300
  UNKNOWN_Z,  // 9_8 unknown 300
  UNKNOWN_Z,  // 9_9 unknown 300
  UNKNOWN_Z,  // 9_10 unknown 300
  UNKNOWN_Z,  // 9_11 unknown 300
  UNKNOWN_Z,  // 9_12 unknown 300
  UNKNOWN_Z,  // 9_13 unknown 300
  UNKNOWN_Z,  // 9_14 unknown 300
  UNKNOWN_Z,  // 9_15 unknown 300
  UNKNOWN_Z,  // 9_16 unknown 300
  UNKNOWN_Z,  // 9_17 unknown 300
  UNKNOWN_Z,  // 9_18 unknown 300
  UNKNOWN_Z,  // 9_19 unknown 300
  UNKNOWN_Z,  // 9_20 unknown 300
  UNKNOWN_Z,  // 9_21 unknown 300
  UNKNOWN_Z,  // 9_22 unknown 300
  UNKNOWN_Z,  // 9_23 unknown 300
  UNKNOWN_Z,  // 9_24 unknown 300
  UNKNOWN_Z,  // 9_25 unknown 300
  UNKNOWN_Z,  // 9_26 unknown 300
  UNKNOWN_Z,  // 9_27 unknown 300
  UNKNOWN_Z,  // 9_28 unknown 300
  UNKNOWN_Z,  // 9_29 unknown 300
  UNKNOWN_Z,  // 9_30 unknown 300
  UNKNOWN_Z,  // 9_31 unknown 300
  UNKNOWN_Z,  // 9_32 unknown 300
  UNKNOWN_Z,  // 9_33 unknown 300
  UNKNOWN_Z,  // 9_34 unknown 300
  UNKNOWN_Z,  // 9_35 unknown 300
  UNKNOWN_Z,  // 9_36 unknown 300
  UNKNOWN_Z,  // 9_37 unknown 300
  UNKNOWN_Z,  // 9_38 unknown 300
  UNKNOWN_Z,  // 9_39 unknown 300
  UNKNOWN_Z,  // 9_40 unknown 300
  UNKNOWN_Z,  // 9_41 unknown 300
  UNKNOWN_Z,  // 9_42 unknown 300
  UNKNOWN_Z,  // 9_43 unknown 300
  UNKNOWN_Z,  // 9_44 unknown 300
  UNKNOWN_Z,  // 9_45 unknown 300
  UNKNOWN_Z,  // 9_46 unknown 300
  UNKNOWN_Z,  // 9_47 unknown 300
  UNKNOWN_Z,  // 9_48 unknown 300
  UNKNOWN_Z,  // 9_49 unknown 300
  UNKNOWN_Z,  // 6c_1 unknown 300
  UNKNOWN_Z,  // 6c_2 unknown 300
  UNKNOWN_Z,  // 7c_1 unknown 300
  UNKNOWN_Z,  // 8c_1 unknown 300
  UNKNOWN_Z,  // 8c_2 unknown 300
  UNKNOWN_Z,  // 8c_3 unknown 300
  UNKNOWN_Z,  // 8c_4 unknown 300
  UNKNOWN_Z}; // 8c_5 unknown 300

char knotname [NUM_KNOTS][5] = {
  "0_1", // 0
  "3_1", // 1
  "4_1", // 2
  "5_1", // 3
  "5_2", // 4
  "6_1", // 5
  "6_2", // 6
  "6_3", // 7
  "7_1", // 8
  "7_2", // 9
  "7_3", // 10
  "7_4", // 11
  "7_5", // 12
  "7_6", // 13
  "7_7", // 14
  "8_1", // 15
  "8_2", // 16
  "8_3", // 17
  "8_4", // 18
  "8_5", // 19
  "8_6", // 20
  "8_7", // 21
  "8_8", // 22
  "8_9", // 23
  "8_10", // 24
  "8_11", // 25
  "8_12", // 26
  "8_13", // 27
  "8_14", // 28
  "8_15", // 29
  "8_16", // 30
  "8_17", // 31
  "8_18", // 32
  "8_19", // 33
  "8_20", // 34
  "8_21", // 35
  "9_1", // 36
  "9_2", // 37
  "9_3", // 38
  "9_4", // 39
  "9_5", // 40
  "9_6", // 41
  "9_7", // 42
  "9_8", // 43
  "9_9", // 44
  "9_10", // 45
  "9_11", // 46
  "9_12", // 47
  "9_13", // 48
  "9_14", // 49
  "9_15", // 50
  "9_16", // 51
  "9_17", // 52
  "9_18", // 53
  "9_19", // 54
  "9_20", // 55
  "9_21", // 56
  "9_22", // 57
  "9_23", // 58
  "9_24", // 59
  "9_25", // 60
  "9_26", // 61
  "9_27", // 62
  "9_28", // 63
  "9_29", // 64
  "9_30", // 65
  "9_31", // 66
  "9_32", // 67
  "9_33", // 68
  "9_34", // 69
  "9_35", // 70
  "9_36", // 71
  "9_37", // 72
  "9_38", // 73
  "9_39", // 74
  "9_40", // 75
  "9_41", // 76
  "9_42", // 77
  "9_43", // 78
  "9_44", // 79
  "9_45", // 80
  "9_46", // 81
  "9_47", // 82
  "9_48", // 83
  "9_49", // 84
  "6c_1", // 85
  "6c_2", // 86
  "7c_1", // 87
  "8c_1", // 88
  "8c_2", // 89
  "8c_3", // 90
  "8c_4", // 91
  "8c_5"}; // 92

// from bfacf.cpp

void delete_Edge(CubicLatticeKnotPtr knot, EdgePtr ep);

// directions are specified by an integer from 0 to 5 inclusive
// MOVE_NORTH is 0, MOVE_EAST is 1, etc.  (see bfacf.h for definitions)

// given a direction DIR, opposite direction is opposite [DIR]
int opposite [6] = {MOVE_SOUTH, MOVE_WEST, MOVE_EAST, MOVE_NORTH, MOVE_DOWN, MOVE_UP};

ivector increment_NEWSUD [6] = {// increments associated with various directions
   { 0, 1, 0}, // North
   { 1, 0, 0}, // East
   {-1, 0, 0}, // West
   { 0, -1, 0}, // South
   { 0, 0, 1}, // Up
   { 0, 0, -1}
}; // Down

int turn [6][4] = {// possible directions for +2 moves, 
   // given direction of edge that is being moved 
   {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // North
   {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // East
   {MOVE_NORTH, MOVE_UP, MOVE_SOUTH, MOVE_DOWN}, // West
   {MOVE_WEST, MOVE_UP, MOVE_EAST, MOVE_DOWN}, // South
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}, // Up
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH}
}; // Down


#define BAD  -1

int reflect [6][6] = {// first index is reflection direction, second index is direction to be reflected
   {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // North
   {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // East
   {MOVE_NORTH, MOVE_WEST, MOVE_EAST, MOVE_SOUTH, MOVE_UP, MOVE_DOWN}, // West
   {MOVE_SOUTH, MOVE_EAST, MOVE_WEST, MOVE_NORTH, MOVE_UP, MOVE_DOWN}, // South
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}, // Up
   {MOVE_NORTH, MOVE_EAST, MOVE_WEST, MOVE_SOUTH, MOVE_DOWN, MOVE_UP}
}; // Down

int kross [6][6] = {
   {BAD, MOVE_DOWN, MOVE_UP, BAD, MOVE_EAST, MOVE_WEST}, // North
   {MOVE_UP, BAD, BAD, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
   {MOVE_DOWN, BAD, BAD, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
   {BAD, MOVE_UP, MOVE_DOWN, BAD, MOVE_WEST, MOVE_EAST}, // South
   {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, BAD, BAD}, // Up
   {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, BAD, BAD}
}; // Down

int rotate90 [6][6] = {// first index is direction of axis of rotation, second index is direction to be rotated
   {MOVE_NORTH, MOVE_DOWN, MOVE_UP, MOVE_SOUTH, MOVE_EAST, MOVE_WEST}, // North
   {MOVE_UP, MOVE_EAST, MOVE_WEST, MOVE_DOWN, MOVE_SOUTH, MOVE_NORTH}, // East
   {MOVE_DOWN, MOVE_EAST, MOVE_WEST, MOVE_UP, MOVE_NORTH, MOVE_SOUTH}, // West
   {MOVE_NORTH, MOVE_UP, MOVE_DOWN, MOVE_SOUTH, MOVE_WEST, MOVE_EAST}, // South
   {MOVE_WEST, MOVE_NORTH, MOVE_SOUTH, MOVE_EAST, MOVE_UP, MOVE_DOWN}, // Up
   {MOVE_EAST, MOVE_SOUTH, MOVE_NORTH, MOVE_WEST, MOVE_UP, MOVE_DOWN}
}; // Down

bool perp [6][6] = {// given direction DIRA and DIRB, perp [DIRA][DIRB] is true 
   // iff DIRA perpendicular to DIRB
   {false, true, true, false, true, true}, // North
   {true, false, false, true, true, true}, // East
   {true, false, false, true, true, true}, // West
   {false, true, true, false, true, true}, // South
   {true, true, true, true, false, false}, // Up
   {true, true, true, true, false, false}
}; // Down

bool anti [6][6] = {// given direction DIRA and DIRB, anti [DIRA][DIRB] is true 
   // iff DIRA antiparallel to DIRB
   {false, false, false, true, false, false}, // North
   {false, false, true, false, false, false}, // East
   {false, true, false, false, false, false}, // West
   {true, false, false, false, false, false}, // South
   {false, false, false, false, false, true}, // Up
   {false, false, false, false, true, false}
}; // Down

// names of the directions
char clk_dir_name [6][6] = {"North", "East", "West", "South", "Up", "Down"};

bool add_edge_to_knot(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, ivector start, ivector end)
{
   EdgePtr ep = (EdgePtr) calloc(1, sizeof (Edge)); // create a new edge

   sub_ivector(ep->increment, end, start); // increment is (end - start)
   ep->dir = clk_direction(ep->increment); // direction associated with increment
   copy_ivector(ep->start, start); // starting location of this edge
   ep->frozen = false; // edge is free to move
   ep->comp = (void *) comp;

   if (comp->first_edge)
   { // some edges have already been added
      comp->last_edge->next = ep;
      ep->prev = comp->last_edge;
      comp->last_edge = ep;
   }
   else
   { // this is the first edge
      comp->first_edge = comp->last_edge = ep;
   }

   comp->first_edge->prev = comp->last_edge;
   comp->last_edge->next = comp->first_edge;
   comp->nedges++;
   clkp->nedges_total++;

   return clk_check_increment(ep->increment);
}

CubicLatticeKnotPtr bfacf_input_start_configuration(FILE *fp)
{
   return bfacf_input_start_configuration(fp, false, DEFAULT_POOLSIZE);
}

CubicLatticeKnotPtr bfacf_input_start_configuration(FILE *fp, bool ignore_mid)
{
   return bfacf_input_start_configuration(fp, ignore_mid, DEFAULT_POOLSIZE);
}

CubicLatticeKnotPtr bfacf_input_start_configuration(FILE *fp, bool ignore_mid, int poolsize)
{
   char in_line [102];
   int line_number = 1;
   int number_vertices_read = 0;
   ivector min, max;

   CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr) calloc(1, sizeof (CubicLatticeKnot));
   ivector first_vertex, vertex, last_vertex;
   bool start_new_comp = true;
   ComponentCLKPtr comp;
   static int comp_ID = 0;

   while (fgets(in_line, 100, fp))
   {
      // a valid start configuration is 3 integers per line
      int num_values_read = sscanf(in_line, "%d%d%d", vertex, vertex + 1, vertex + 2);
      if (num_values_read < 3)
      {
         start_new_comp = true;
         if (comp) // close off the previous component
            add_edge_to_knot(clkp, comp, last_vertex, first_vertex);

         continue;
      }
      else if (num_values_read != 3)
      {
         fprintf(stderr, "error on line #%d found: %d >%s<\n", line_number, num_values_read, in_line);
         return (CubicLatticeKnotPtr) NULL;
      }

      if (start_new_comp)
      { // read first vertex
         copy_ivector(first_vertex, vertex); // remember the first vertex for later
         start_new_comp = false;
         comp = (ComponentCLKPtr) calloc(1, sizeof (ComponentCLK));
         comp->clkp = clkp;
         comp->minedges = clk_minedges;
         comp->maxedges = clk_maxedges;
         //comp->ID = comp_ID++;  now done later in  bfacf_init_pool_and_lattice()
         if (clkp->fcomp)
         { // clkp already has some components
            clkp->lcomp->next = comp;
            comp->prev = clkp->lcomp;
            clkp->lcomp = comp;
            clkp->ncomps++;
         }
         else
         {
            clkp->fcomp = clkp->lcomp = comp;
            clkp->ncomps = 1;
         }

         if (clkp->ncomps == 1)
         {
            copy_ivector(min, vertex);
            copy_ivector(max, vertex);
         }
      }
      else
      {
         max_ivector(max, max, vertex);
         min_ivector(min, min, vertex);
         if (!add_edge_to_knot(clkp, comp, last_vertex, vertex))
         {
            fprintf(stderr, "bad increment of (%d, %d, %d) seen on line #%d\n",
                    comp->last_edge->increment [0],
                    comp->last_edge->increment [1],
                    comp->last_edge->increment [2],
                    line_number);
            return (CubicLatticeKnotPtr) NULL;
         }
      }
      copy_ivector(last_vertex, vertex);
      ++line_number;
   }
   add_edge_to_knot(clkp, comp, vertex, first_vertex);
   /*
   if (number_vertices_read < 4) {
     fprintf (stderr, " *** Knot has too few (%d) vertices\n", number_vertices_read);
     return (CubicLatticeKnotPtr) NULL;
   }
    */

   ivector range;
   sub_ivector(range, max, min);
   if (range [0] > LATTICE_SIZE - 3 || range [1] > LATTICE_SIZE - 3 || range [2] > LATTICE_SIZE - 3)
   {
      fprintf(stderr, " *** Knot has excessive extent in one or more dimensions\n");
      return (CubicLatticeKnotPtr) NULL;
   }

   if (poolsize < (108 * clkp->nedges_total) / 100 + 108)
   {
      fprintf(stderr, "poolsize too small\n");
      free(clkp);
      return (CubicLatticeKnotPtr) NULL;
   }

   bfacf_init_pool_and_lattice(clkp, ignore_mid, poolsize, min, max);

   return clkp;
}

static char *clk_recycled_lattice = (char *) NULL;

void clk_set_free_recycled_lattice(char *lat) { }

void bfacf_init_pool_and_lattice(CubicLatticeKnotPtr clkp, bool ignore_mid, int poolsize, ivector min, ivector max)
{
   ComponentCLKPtr comp = clkp->fcomp;

   clkp->auto_recentre = true;

   // create a `swimming space' (a pool) for our edges
   // this is to avoid frequent memory allocation / deallocation

   clkp->poolsize = poolsize;
   clkp->edgepool = (EdgePtr *) calloc(clkp->poolsize, sizeof (EdgePtr));
   int index = 0;
   int ID = 0;

   comp = clkp->fcomp;
   while (comp)
   {
      comp->ID = ID++;
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (index > poolsize)
         {
            printf("yikes! %d %d\n", index, poolsize);
            exit(11);
         }
         clkp->edgepool [index++] = ep;
         ep = ep->next;
      }
      comp = comp->next;
   }

   while (index < poolsize)
      clkp->edgepool [index++] = (EdgePtr) calloc(1, sizeof (Edge)); // new, unused edge

   for (index = 0; index < clkp->poolsize; index++)
   {
      clkp->edgepool [index]->locpool = index;
      clkp->edgepool [index]->ID = index;
   }

   clkp->lattice = (char *) calloc(LATTICE_TOTAL_SIZE, sizeof (char));
   ivector mid;
   midpoint(mid, min, max);

   if (ignore_mid)
      set_ivector(mid, 0, 0, 0);

   init_lattice(clkp, mid);


   comp = clkp->fcomp;
   while (comp)
   {
      ivector closing_incr;
      sub_ivector(closing_incr, comp->first_edge->start, comp->last_edge->start);
      // if the edge that would close the knot is not valid,
      // then we have an open-ended string.  
      if (!clk_check_increment(closing_incr))
      {
         comp->flags |= COMPONENT_CLK_FLAG_OPEN;
         delete_Edge(clkp, comp, comp->last_edge);
         freezeEdge(clkp, comp->first_edge);
         freezeEdge(clkp, comp->last_edge);
         comp->first_edge->frozen = comp->last_edge->frozen = true;
      }

      comp = comp->next;
   }
}

bool enlarge_pool(CubicLatticeKnotPtr knot)
{
   return false;
   //  fprintf (stderr, "enlarge_pool() not implemented yet! Knot has %d edges\n", knot->nedges);
   //  exit (102);
}

bool reject_via_energy(double previous_energy, double current_energy)
{
   if (current_energy > previous_energy)
   {
      double p = exp(previous_energy - current_energy); // probability of acceptance
      if (rand_uniform() > p)
      {
         return true; // reject it
      }
   }

   return false;
}

// the three moves below assume that that move is legal locally, but maybe not globally

bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep, ivector increment, int dir)
{
   // increment is the direction of the move
   // dir is the direction associated with this increment

   if (knot->nedges_total > knot->poolsize - 3)
   {
      enlarge_pool(knot);
      return false;
   }
   if (comp->nedges > comp->maxedges) return false;

   // location that edge might want to move to must be empty

   ivector test_locationA, test_locationB;
   add_ivector(test_locationA, ep->start, increment);

   if (clk_check_for_edge_hits(knot, dir, test_locationA))
   {
      if (knot->auto_recentre)
      {
         if (recentre_knot_in_lattice(knot)) // returns true if knot can't be recentered
            return false;
         add_ivector(test_locationA, ep->start, increment); // need to recalculate this
      }
      else
         return false;
   }

   if (knot->lattice [lat(test_locationA)] == OCCUPIED) return false;
   add_ivector(test_locationB, ep->next->start, increment);
   if (knot->lattice [lat(test_locationB)] == OCCUPIED) return false;

   // rejection by energy should go here (+2 move)
   if (knot->filter_by_energy)
   {
      double E_before =
              energy(ep->start, knot, ep->next->next, ep->prev->prev) +
              energy(ep->next->start, knot, ep->next->next->next, ep->prev);
      double E_after =
              energy(ep->start, test_locationB) +
              energy(ep->start, knot, ep->next, ep->prev->prev) +
              energy(test_locationA, knot, ep->next, ep->prev) +
              energy(test_locationB, knot, ep->next->next, ep->prev) +
              energy(ep->next->start, knot, ep->next->next->next, ep->prev);
      if (reject_via_energy(E_before, E_after)) return false;
   }

   // choose the first two unused edges in edge pool
   EdgePtr ep_prev = knot->edgepool [knot->nedges_total]; // note: not nedges - 1 as in above
   EdgePtr ep_next = knot->edgepool [knot->nedges_total + 1]; // note: not nedges - 1 as in above
   ep_prev->comp = ep_next->comp = (void *) comp;

   // adjust adjacency pointers
   ep_prev->prev = ep->prev;
   ep_prev->next = ep;
   ep_next->prev = ep;
   ep_next->next = ep->next;
   ep->prev->next = ep_prev;
   ep->next->prev = ep_next;
   ep->prev = ep_prev;
   ep->next = ep_next;

   copy_ivector(ep_prev->increment, increment);
   negate_ivector(ep_next->increment, increment);
   ep_prev->dir = dir;
   ep_next->dir = opposite [dir];


   // increase by 2 number of edges in comp and in knot
   comp->nedges += 2;
   knot->nedges_total += 2;

   // update positions 

   EdgePtr ep2 = ep_prev;
   for (int i = 0; i < 3; i++)
   {
      add_ivector(ep2->start, ep2->prev->start, ep2->prev->increment);
      ep2 = ep2->next;
   }

   // update lattice 

   knot->lattice [lat(ep->start)] = OCCUPIED;
   knot->lattice [lat(ep->next->start)] = OCCUPIED;


   knot->success_plus2++;
#ifdef TESTING
   knot->p2dir [dir]++;
#endif
   return true;
}

bool perform_plus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // this function assumes that a move in any of the four directions
   // is possible locally (ignoring global self-intersections)

   // choose a direction at random
   int dir = turn [ep->dir][rand_integer(0, 4)];
   return perform_plus2_move(knot, comp, ep, increment_NEWSUD [dir], dir);
}

bool perform_plus2_move_alt(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   int dir;
   // this function is similar to perform_plus2_move() above, but
   // assumes that a move in at most three of the four directions
   // is possible locally (ignoring global self-intersections)

   // in order for the move to be possible locally, direction of
   // the move cannot be in the same direction as the next edge,
   // and not in the opposite direction as the previous edge

   do
   {
      dir = turn [ep->dir][rand_integer(0, 4)];
   }
   while (dir == ep->next->dir || dir == opposite [ep->prev->dir]);

   return perform_plus2_move(knot, comp, ep, increment_NEWSUD [dir], dir);
}

bool perform_0_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // type 0 moves are only possible if this edge and the next edge
   // are perpendicular

#ifdef TESTING
   if (!perp [ep->dir][ep->next->dir])
   {
      fprintf(stderr,
              "\nperform_0_move (): edges are not perpendicular as expected!  Panicking!!\n\n");
      fprintf(stderr, "offending edge is %d (%d <--> %d <--> %d)\n",
              ep->ID, ep->prev->ID, ep->ID, ep->next->ID);
      for (int i = 0; i < 22; i++)
         fprintf(stderr, "%d ", knot->edgepool [i]->ID);
      exit(102);
   }
#endif

   if (ep->frozen || ep->next->frozen) return false; // NEW!!

   // also the new lattice location needs to be empty
   ivector test_location;
   add_ivector(test_location, ep->start, ep->next->increment);
   if (knot->lattice [lat(test_location)] == OCCUPIED) return false;

   // rejection by energy should go here (0 move)

   if (knot->filter_by_energy)
   {
      if (reject_via_energy(
              energy(ep->next->start, knot, ep->next->next->next, ep->prev),
              energy(test_location, knot, ep->next->next->next, ep->prev))) return false;
   }

   ivector tmp;
   // simply swap increments between this edge and next
   copy_ivector(tmp, ep->increment);
   copy_ivector(ep->increment, ep->next->increment);
   copy_ivector(ep->next->increment, tmp);

   // clear lattice at old vertex location
   knot->lattice [lat(ep->next->start)] = EMPTY;

   int itmp;
   itmp = ep->dir;
   ep->dir = ep->next->dir;
   ep->next->dir = itmp;

   // update position and lattice
   add_ivector(ep->next->start, ep->start, ep->increment);
   knot->lattice [lat(ep->next->start)] = OCCUPIED;

   knot->success_0++;

   return true;
}

void delete_Edge(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   // adjust adjacent edges' adjacency pointers
   ep->prev->next = ep->next;
   ep->next->prev = ep->prev;

   // swap this edge with the last edge in the edgepool
   EdgePtr tmp = knot->edgepool [knot->nedges_total - 1];
   tmp->locpool = ep->locpool;
   knot->edgepool [ep->locpool] = tmp;
   knot->edgepool [knot->nedges_total - 1] = ep;
   ep->locpool = knot->nedges_total - 1;

   // reduce by 1 the number of edges in knot and in comp
   --knot->nedges_total;
   --comp->nedges;

#ifdef TESTING
   if (comp->nedges < 4)
   {
      fprintf(stderr, "Absurd number of edges (%d) found in knot!  Panicking!\n", comp->nedges);
      fflush(stderr);
      exit(1);
   }
#endif

   // update info contained in `comp', in the event we've 
   // deleted first or last edge in that component

   if (comp->first_edge == ep)
      comp->first_edge = ep->next;
   else if (comp->last_edge == ep)
      comp->last_edge = ep->prev;
}

bool perform_minus2_move(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   if (comp->nedges < 5) return false;
   if (comp->nedges < comp->minedges) return false;
   if (ep->next->frozen || ep->prev->frozen) return false;

#ifdef TESTING

   // type -2 moves are only possible if previous edge and next edge are anti-parallel
   if (!anti [ep->prev->dir][ep->next->dir])
   {
      fprintf(stderr,
              "\nAdjacent edges are not anti-parallel as expected!  Panicking!!\n\n");
      exit(102);
   }

   // also edge and next edge need to be perpendicular
   if (!perp [ep->dir][ep->next->dir])
   {
      fprintf(stderr, "\nEdges are not perpendicular as expected!  Panicking!!\n\n");
      exit(102);
   }

   knot->pm2dir [ep->next->dir]++;
#endif

   // rejection by energy should go here (-2 move)
   if (knot->filter_by_energy)
   {
      double E_before =
              energy(ep->prev->start, knot, ep->next, ep->prev->prev->prev) +
              energy(ep->start, knot, ep->next->next, ep->prev->prev) +
              energy(ep->next->start, knot, ep->next->next->next, ep->prev->prev) +
              energy(ep->next->next->start, knot, ep->next->next->next->next, ep->prev->prev);
      double E_after =
              energy(ep->prev->start, knot, ep->next->next->next, ep->prev->prev->prev) +
              energy(ep->next->next->start, knot, ep->next->next->next->next, ep->prev->prev);
      if (reject_via_energy(E_before, E_after)) return false;
   }

   // clear lattice 
   knot->lattice [lat(ep->start)] = EMPTY;
   knot->lattice [lat(ep->next->start)] = EMPTY;

   // for this case, no increments need to be changed or specified
   delete_Edge(knot, comp, ep->prev);
   delete_Edge(knot, comp, ep->next);

   // update position
   add_ivector(ep->start, ep->prev->start, ep->prev->increment);

   knot->success_minus2++;
   return true;
}

bool bad_pivot_choice(int &dir, EdgePtr ep1, EdgePtr ep2)
{
   if (ep1 == ep2) return true;
   // if (ep1->comp != ep2->comp) return true;
   // Multi-component links are simply banned outright before
   // we get to this stage. The components quickly separate 
   // and drift away from each other.
   ivector diff;
   sub_ivector(diff, ep1->start, ep2->start);
   int m = MAX(abs(diff [0]), abs(diff [1]));
   m = MAX(m, abs(diff [2]));
   //printf ("%d %d %d\n", diff [0], diff [1], diff [2]); fflush (stdout);
   if (m < 4) return true;

   int index;
   m = 0;
   for (int i = 0; i < 3; i++)
   {
      if (diff [i] == 0)
         m++;
      else
         index = i;
   }
   if (m != 2) return true;
   switch (index)
   {
      case 0:
         if (diff [0] < 0)
            dir = MOVE_WEST;
         else
            dir = MOVE_EAST;
         break;
      case 1:
         if (diff [1] < 0)
            dir = MOVE_SOUTH;
         else
            dir = MOVE_NORTH;
         break;
      case 2:
         if (diff [2] < 0)
            dir = MOVE_DOWN;
         else
            dir = MOVE_UP;
         break;
   }
   return false;
}

void set_sub_path(CubicLatticeKnotPtr knot, EdgePtr ep1, EdgePtr ep2, int value)
{ // NOTE: sets start of both ep1 and ep2
   EdgePtr ep = ep1;
   while (ep != ep2->next)
   {
      knot->lattice [lat(ep->start)] = value;
      ep = ep->next;
   }
}

int rtemp [6];

bool perform_move_pivot(CubicLatticeKnotPtr knot)
{
   static unsigned int call = 0;
   double p = rand_uniform(); // uniform number between 0 and 1
   bool value = false; // assume move fails
   //clk_allocation_alt_lattice (knot);
   //clk_fill_path_alt (knot, true);
   int index;
   EdgePtr ep1, ep2;
   // choose two edges uniformly at random from the edge pool
   int kount = 0;
   do
   {
      ep1 = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
      ep2 = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
      kount++;
   }
   while (bad_pivot_choice(index, ep1, ep2));
   //printf ("%d ", kount);fflush (stdout);
   EdgePtr race1 = ep1;
   EdgePtr race2 = ep2;
   while (!(race1 != ep2 || race2 != ep1))
   {
      race1 = race1->next;
      race2 = race2->next;
   }
   if (race2 == ep1)
   {
      ep1 = ep2;
      ep2 = race2;
      index = opposite [index];
   }

   int axis_dir = index;
   if (axis_dir < 0)
   {
      fprintf(stderr, "\n\n\nbad axis!  %d  (%d %d %d) --- (%d %d %d)\n",
              axis_dir,
              ep1->start [0], ep1->start [1], ep1->start [2],
              ep2->start [0], ep2->start [1], ep2->start [2]);
      exit(333);
   }
   int numb90 = rand_integer(1, 8); // number of 90 degree rotations, either 1, 2 or 3 or reflection direction 4, 5, 6 or 7

   if (numb90 < 4)
   {
      for (int i = 0; i < 6; i++)
      {
         int k = numb90;
         rtemp [i] = i;
         while (k--)
            rtemp [i] = rotate90 [axis_dir][rtemp [i]];
      }
   }
   else
   {
      int reflection_direction = turn [axis_dir][numb90 - 4];
      for (int i = 0; i < 6; i++)
         rtemp [i] = reflect [reflection_direction][i];
   }
   // see if we can get from ep1 to ep2
   EdgePtr ep = ep1;
   ivector loc;
   ivector par1;
   copy_ivector(par1, ep2->start); // paranoia check
   copy_ivector(loc, ep1->start);
   set_sub_path(knot, ep1, ep2, EMPTY); // clear current path
   while (ep != ep2)
   {
      int dir = rtemp [ep->dir];
      add_ivector(loc, loc, increment_NEWSUD [dir]);
      ep->scratch = dir;
      if (knot->lattice [lat(loc)] != EMPTY)
      {
         set_sub_path(knot, ep1, ep2, OCCUPIED); // refill path
         return false;
      }
      ep = ep->next;
   }

   // if we got here, path is clear
   // update directions, increments and start locations

   copy_ivector(par1, ep2->start); // paranoia check
   ep = ep1;
   knot->lattice [lat(ep->start)] = OCCUPIED;
   while (ep != ep2)
   {
      ep->dir = ep->scratch;
      copy_ivector(ep->increment, increment_NEWSUD [ep->dir]);
      add_ivector(ep->next->start, ep->start, ep->increment);
      knot->lattice [lat(ep->next->start)] = OCCUPIED;
      ep = ep->next;
   }

#ifdef DEBUG_PIVOT
   static int pivot_debug = 0;
   clk_validate(knot, "perform_move_pivot");
   if (!equal_ivector(par1, ep2->start))
   {
      fprintf(stderr, "pivot panic!!! rotating about %d %s by %d 90 degree turns\n\n",
              axis_dir, clk_dir_name [axis_dir], numb90);
      for (int i = 0; i < 6; i++)
         fprintf(stderr, "%s becomes %s (%d)\n", clk_dir_name [i], clk_dir_name [rtemp [i]], rtemp [i]);
      panic_exit("pivot panic\n\n");
   }
   //set_sub_path (knot, ep1, ep2, OCCUPIED);  // refill path
   if (++pivot_debug % 20000 == 0)
   {
      printf("%d ", pivot_debug / 20000);
      fflush(stdout);
   }
#endif

   if (call++ % 10 == 0) recentre_knot_in_lattice(knot);
   return true;
}

bool perform_move(CubicLatticeKnotPtr knot)
{
   if (knot->nfrozen == knot->nedges_total) return false;

   // choose an edge uniformly at random from the edge pool

   EdgePtr ep = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;

#ifdef CHECK_CONFIG_not_used
   ivector eploc;
   copy_ivector(eploc, ep->start);
#endif

   // The types of moves possible depend on the chosen edge, ep, and
   // the two adjacent edges, ep->prev and ep->next.

   double p = rand_uniform(); // uniform number between 0 and 1
//   printf("%f ", p);
   bool value = false; // assume move fails

   // NOTE: reordering the cases below may increase performance slightly
   // but Case 1 needs to be tested before Case 3.

   // Case 1: chosen edge is parallel to both adjacent edges,
   // all four moves are +2 moves.

   // directions of adjacent edges same as chosen edge

   if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
   {
      if (p < comp->p_4p2)
         value = perform_plus2_move(knot, comp, ep);
   }

      // Case 2: chosen edge is perpendicular to both adjacent edges
      // and those adjacent edges are anti-parallel to each other. 
      // One move is a -2 move, other three are +2 moves.

   else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
   {
      if (p < comp->p_minus2)
         value = perform_minus2_move(knot, comp, ep);
      else if (p < comp->p_m23p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

      // Case 3: chosen edge is perpendicular to one adjacent edge
      // and parallel to the other adjacent edge.
      // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
      // because the case of both being parallel has been ruled out above.
      // One move is a 0 move, other three are +2 moves.

      // NOTE: test this before Case 1 above

   else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
   {
      if (p < comp->p_0)
      {
         if (ep->dir == ep->prev->dir)
            value = perform_0_move(knot, comp, ep);
         else
            value = perform_0_move(knot, comp, ep->prev);
      }
      else if (p < comp->p_03p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

      // Case 4: Only possibility remaining: chosen edge is perpendicular
      // to both adjacent edges, and the adjacent edges are not anti-parallel.
      // Two moves are 0 moves and two moves are +2 moves.

   else
   {
      if (p < comp->p_2p0)
      {
         if (p < comp->p_0)
            value = perform_0_move(knot, comp, ep->prev);
         else
            value = perform_0_move(knot, comp, ep);
      }
      else if (p < comp->p_2p02p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

#ifdef CHECK_CONFIG_NOT_USED
   extern void bfacf_check_config(ivector);
   bfacf_check_config(eploc);
#endif

#ifdef TESTING2    // Warning: this changes the algorithm from roughly constant time
   //          to roughly linear time (in number of edges)
   if (value)
      clk_check_increments(knot);

#endif

   return value;
}

bool perform_crankshaft(CubicLatticeKnotPtr knot, ComponentCLKPtr comp, EdgePtr ep)
{
   if (comp->nedges < 5) return false;
   int dir = ep->prev->dir; // remember the direction of the previous edge
   perform_minus2_move(knot, comp, ep); // get rid of the edge
   if (perform_plus2_move(knot, comp, ep))
      return true;
   else // the above has failed, move the edge back to where it was
      return perform_plus2_move(knot, comp, ep, increment_NEWSUD [dir], dir);
}

bool perform_move_kjc(CubicLatticeKnotPtr knot)
{
   // choose an edge uniformly at random from the edge pool
   if (knot->nfrozen == knot->nedges_total) return false;
   EdgePtr ep = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;

   // The types of moves possible depend on the chosen edge, ep, and
   // the two adjacent edges, ep->prev and ep->next.

   double p = rand_uniform(); // uniform number between 0 and 1
   bool value = false; // assume move fails

   if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
   {
      value = perform_crankshaft(knot, comp, ep);
   }
   else if (perp [ep->dir][ep->prev->dir] || perp [ep->dir][ep->next->dir])
   {
      if (rand_uniform() < 0.5)
      {
         if (perp [ep->dir][ep->prev->dir])
            value = perform_0_move(knot, comp, ep->prev);
      }
      else
      {
         if (perp [ep->dir][ep->next->dir])
            value = perform_0_move(knot, comp, ep);
      }
   }

   return value;
}

bool perform_move_move0_only(CubicLatticeKnotPtr knot)
{
   // choose an edge uniformly at random from the edge pool

   EdgePtr ep = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;

   // this might be somewhat faster than using BFACF and setting probabilities to 0 1 0

   double p = rand_uniform(); // uniform number between 0 and 1
   bool value = false; // assume move fails

   if (perp [ep->dir][ep->prev->dir] || perp [ep->dir][ep->next->dir])
   {
      if (rand_uniform() < 0.5)
      {
         if (perp [ep->dir][ep->prev->dir])
            value = perform_0_move(knot, comp, ep->prev);
      }
      else
      {
         if (perp [ep->dir][ep->next->dir])
            value = perform_0_move(knot, comp, ep);
      }
   }

   return value;
}

bool type0_moves_possible(CubicLatticeKnotPtr clkp)
{
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (perp [ep->dir][ep->next->dir])
         {
            ivector test_location;
            add_ivector(test_location, ep->start, ep->next->increment);
            if (clkp->lattice [lat(test_location)] == EMPTY) return true;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

int number_type0_moves(CubicLatticeKnotPtr clkp)
{
   int number;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      int number = 0;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (perp [ep->dir][ep->next->dir])
         {
            ivector test_location;
            add_ivector(test_location, ep->start, ep->next->increment);
            if (clkp->lattice [lat(test_location)] == EMPTY) ++number;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return number;
}

bool type_minus2_moves_possible(CubicLatticeKnotPtr clkp)
{
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (perp [ep->dir][ep->next->dir] && anti [ep->dir][ep->next->next->dir]) return true;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

int number_type_minus2_moves(CubicLatticeKnotPtr clkp)
{
   ComponentCLKPtr comp = clkp->fcomp;
   int number = 0;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (perp [ep->dir][ep->next->dir] && anti [ep->dir][ep->next->next->dir]) ++number;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return number;
}

int count_plus2_moves(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   // in order for the move to be possible locally, direction of
   // the move cannot be in the same direction as the next edge,
   // and not in the opposite direction as the previous edge
   int kount = 0;

   for (int i = 0; i < 4; i++)
   {
      int dir = turn [ep->dir][i];
      if (dir != ep->next->dir && dir != opposite [ep->prev->dir])
      {
         // move is possible locally, but maybe not globally
         ivector test_location1, test_location2;
         add_ivector(test_location1, ep->start, increment_NEWSUD [dir]);
         add_ivector(test_location2, ep->next->start, increment_NEWSUD [dir]);
         // move is possible globally if both test locations are empty
         if (clkp->lattice [lat(test_location1)] == EMPTY && clkp->lattice [lat(test_location2)] == EMPTY) ++kount;
      }
   }

   return kount;
}

void bfacf_number_moves(CubicLatticeKnotPtr clkp, int &nm2, int &n0, int &np2)
{
   nm2 = n0 = np2 = 0;
   if (clkp->nedges_total == 4)
   {
      np2 = 12;
      return;
   }
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {

         np2 += count_plus2_moves(clkp, ep);

         if (perp [ep->dir][ep->next->dir])
         {
            ivector test_location;
            add_ivector(test_location, ep->start, ep->next->increment);
            if (clkp->lattice [lat(test_location)] == EMPTY) ++n0;
            if (anti [ep->dir][ep->next->next->dir]) ++nm2;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
}

int bfacf_number_sticks(CubicLatticeKnotPtr clkp)
{
   int dirchanges = 0;
   if (!clkp) return 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (ep->dir != ep->next->dir)
            dirchanges++;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return dirchanges;
}

EdgePtr clk_get_edge(CubicLatticeKnotPtr clkp, ivector start)
{
   ComponentCLKPtr comp = clkp->fcomp;
   add_ivector(start, start, clkp->loffset);
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(ep->start, start))
         {
            return ep;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL;
}

void swapEdgesPool(CubicLatticeKnotPtr clkp, int loc1, int loc2)
{
   if (loc1 == loc2) return;
   EdgePtr ept1 = clkp->edgepool [loc1];
   EdgePtr ept2 = clkp->edgepool [loc2];
   clkp->edgepool [loc1] = ept2;
   clkp->edgepool [loc2] = ept1;
   ept1->locpool = loc2;
   ept2->locpool = loc1;
}

void freezeEdge(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->locpool < clkp->nfrozen) // edge is already frozen
      return;
   swapEdgesPool(clkp, ep->locpool, clkp->nfrozen);
   clkp->nfrozen++;
}

bool freezeEdge(CubicLatticeKnotPtr clkp, ivector start)
{
   // freeze edge in `clkp' that starts at `start'
   if (!clkp) return false;
   EdgePtr ep = clk_get_edge(clkp, start);
   if (!ep) return false;
   freezeEdge(clkp, ep);
   ep->frozen = true;
   return true;
}

void freezeEdge(EdgePtr ep)
{
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr) comp->clkp;
   freezeEdge(clkp, ep);
}

void thawEdge(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->locpool >= clkp->nfrozen) // edge is already thawed
      return;
   swapEdgesPool(clkp, ep->locpool, clkp->nfrozen - 1);
   clkp->nfrozen--;
}

void thawEdge(EdgePtr ep)
{
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   CubicLatticeKnotPtr clkp = (CubicLatticeKnotPtr) comp->clkp;
   thawEdge(clkp, ep);
}

bool freeze_component(CubicLatticeKnotPtr clkp, int ID, bool freeze)
{
   if (!clkp) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      if (ID == comp->ID)
      {
         EdgePtr ep = comp->first_edge;
         for (int i = 0; i < comp->nedges; i++)
         {
            if (freeze)
               freezeEdge(clkp, ep);
            else
               thawEdge(clkp, ep);
            ep = ep->next;
         }
         return true;
      }
      comp = comp->next;
   }
   return false;
}

#ifdef INCLUDED_FROM_KnotPlot

int bfacf_save_stick(int nsticks, CubicLatticeKnotPtr clkp)
{
   int dirchanges = bfacf_number_sticks(clkp);
   if (dirchanges <= nsticks)
   {
      decode_line("bfacf news");
      return dirchanges;
   }
   return 0;
}


#endif

int clk_minedges = 0;
int clk_maxedges = HUGE_NUMBER;

void clk_set_edge_limit(int min, int max)
{ // for newly created components
   clk_minedges = min;
   clk_maxedges = max;
}

void clk_set_nedge_limit(CubicLatticeKnotPtr clkp, int min, int max, int ID)
{
   if (!clkp) return;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      if (ID < 0 || ID == comp->ID)
      {
         comp->minedges = min;
         comp->maxedges = max;
      }
      comp = comp->next;
   }
}

bool perform_move_q(CubicLatticeKnotPtr knot) //cut and paste this code into overridden step function
{
   if (knot->nfrozen == knot->nedges_total) return false;
   if (knot->nedges_total >= knot->q_prob_size) return false; // FIX LATER
   // choose an edge uniformly at random from the edge pool

   EdgePtr ep = knot->edgepool [rand_integer(knot->nfrozen, knot->nedges_total)];
   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;

#ifdef CHECK_CONFIG_not_used
   ivector eploc;
   copy_ivector(eploc, ep->start);
#endif

   // The types of moves possible depend on the chosen edge, ep, and
   // the two adjacent edges, ep->prev and ep->next.

   double p = rand_uniform(); // uniform number between 0 and 1
   bool value = false; // assume move fails

   // NOTE: reordering the cases below may increase performance slightly
   // but Case 1 needs to be tested before Case 3.

   // Case 1: chosen edge is parallel to both adjacent edges,
   // all four moves are +2 moves.

   // directions of adjacent edges same as chosen edge

   if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
   {
      if (p < knot->q_prob [comp->nedges].p_4p2)
         value = perform_plus2_move(knot, comp, ep);
   }

      // Case 2: chosen edge is perpendicular to both adjacent edges
      // and those adjacent edges are anti-parallel to each other. 
      // One move is a -2 move, other three are +2 moves.

   else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
   {
      if (p < knot->q_prob [comp->nedges].p_minus2)
         value = perform_minus2_move(knot, comp, ep);
      else if (p < knot->q_prob [comp->nedges].p_m23p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

      // Case 3: chosen edge is perpendicular to one adjacent edge
      // and parallel to the other adjacent edge.
      // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
      // because the case of both being parallel has been ruled out above.
      // One move is a 0 move, other three are +2 moves.

      // NOTE: test this before Case 1 above

   else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
   {
      if (p < knot->q_prob [comp->nedges].p_0)
      {
         if (ep->dir == ep->prev->dir)
            value = perform_0_move(knot, comp, ep);
         else
            value = perform_0_move(knot, comp, ep->prev);
      }
      else if (p < knot->q_prob [comp->nedges].p_03p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

      // Case 4: Only possibility remaining: chosen edge is perpendicular
      // to both adjacent edges, and the adjacent edges are not anti-parallel.
      // Two moves are 0 moves and two moves are +2 moves.

   else
   {
      if (p < knot->q_prob [comp->nedges].p_2p0)
      {
         if (p < knot->q_prob [comp->nedges].p_0)
            value = perform_0_move(knot, comp, ep->prev);
         else
            value = perform_0_move(knot, comp, ep);
      }
      else if (p < knot->q_prob [comp->nedges].p_2p02p2)
         value = perform_plus2_move_alt(knot, comp, ep);
   }

#ifdef CHECK_CONFIG_NOT_USED
   extern void bfacf_check_config(ivector);
   bfacf_check_config(eploc);
#endif

#ifdef TESTING2    // Warning: this changes the algorithm from roughly constant time
   //          to roughly linear time (in number of edges)
   if (value)
      clk_check_increments(knot);

#endif

   return value;
}

// from clk_util.cpp

bool clk_recombo_limit = false;
int clk_min_arc_length_distance = 0;

int clk_count_edges(EdgePtr ep)
{
   int kount = 1;
   EdgePtr start = ep;
   ep = ep->next;
   while (ep != start)
   {
      kount++;
      ep = ep->next;
   }
   return kount;
}

int clk_validate(CubicLatticeKnotPtr clkp)
{
   for (int i = 0; i < clkp->nedges_total; i++)
      if (clkp->edgepool [i]->locpool != i)
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
         if ((ComponentCLKPtr) ep->comp != comp)
            return 5;
         if (clkp->lattice [lat(ep->start)] != OCCUPIED)
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

void clk_validate(CubicLatticeKnotPtr clkp, char *s)
{
   int code = clk_validate(clkp);
   if (code)
   {
      fprintf(stderr, "fails at %s with code %d\n", s, code);
      exit(102);
   }
   //  else 
   //  fprintf (stderr, "ok at %s\n", s);
}

bool clk_check_increment(ivector incr)
{
   // a valid increment should be +-1 in one direction only
   if (abs(incr [0]) + abs(incr [1]) + abs(incr [2]) == 1)
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
         if (knot->lattice [lat(ep->start)] != OCCUPIED)
         {
            printf("\n\nlattice location %d %d %d for edge with ID %d is not occupied!\n",
                    ep->start [0], ep->start [1], ep->start [2], ep->ID);
            fflush(stdout);
            exit(102);
         }
         //printf ("%d ", ep->ID); fflush (stdout);
         knot->lattice [lat(ep->start)] = EMPTY;

         ep = ep->next;
      }

      if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
      {
         ivector endpoint;
         add_ivector(endpoint, comp->last_edge->start, comp->last_edge->increment);
         knot->lattice [lat(endpoint)] = EMPTY;
      }

      comp = comp->next;
   }

   int ix, iy, iz;
   ivector loc;

   // we've cleared the knot path, so every location in lattice should now be EMPTY

   for (ix = 0; ix < LATTICE_SIZE; ix++)
   { // unnecessesary brace to make VC++ indenter happy
      for (iy = 0; iy < LATTICE_SIZE; iy++)
      {
         for (iz = 0; iz < LATTICE_SIZE; iz++)
         {
            set_ivector(loc, ix, iy, iz);
            if (knot->lattice [lat(loc)] != EMPTY)
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

int clk_arc_length_distance(EdgePtr ep1, EdgePtr ep2)
{
   ComponentCLKPtr comp1 = (ComponentCLKPtr) ep1->comp;
   ComponentCLKPtr comp2 = (ComponentCLKPtr) ep2->comp;
   if (comp1 != comp2)
      return ARC_LENGTH_DISTANCE_INFINITE;

   EdgePtr epn = ep1->next;
   EdgePtr epp = ep1->prev;
   int distance = 1;
   while (epn != ep2 && epp != ep2)
   {
      epn = epn->next;
      epp = epp->prev;
      ++distance;
   }

   return distance;
}

int clk_arc_length_distance(CubicLatticeKnotPtr clkp, ivector s1, ivector s2)
{
   EdgePtr ep1 = clk_get_edge(clkp, s1);
   if (!ep1) return ARC_LENGTH_DISTANCE_NO_EDGE1;
   EdgePtr ep2 = clk_get_edge(clkp, s2);
   if (!ep2) return ARC_LENGTH_DISTANCE_NO_EDGE2;
   return clk_arc_length_distance(ep1, ep2);
}

void clk_get_extent(ivector minbb, ivector maxbb, CubicLatticeKnotPtr knot)
{
   if (!knot) return;
   copy_ivector(minbb, knot->fcomp->first_edge->start);
   copy_ivector(maxbb, minbb);
   ComponentCLKPtr comp = knot->fcomp;

   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         min_ivector(minbb, minbb, ep->start);
         max_ivector(maxbb, maxbb, ep->start);
         ep = ep->next;
      }
      comp = comp->next;
   }
   sub_ivector(minbb, minbb, knot->loffset);
   sub_ivector(maxbb, maxbb, knot->loffset);
}

void clk_bounding_box(int &X, int &Y, int &Z, CubicLatticeKnotPtr knot)
{
   ivector minbb, maxbb;
   clk_get_extent(minbb, maxbb, knot);
   X = maxbb [0] - minbb [0];
   Y = maxbb [1] - minbb [1];
   Z = maxbb [2] - minbb [2];
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
   FILE *fp = (FILE *) NULL;
   //fp = fopen ("out.txt", "w");

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
                    ep->increment [0], ep->increment [1], ep->increment [2]);
            exit(343);
         }
         if (clk_direction(ep->increment) != ep->dir)
         {
            fprintf(stderr, "\n\n clk_check_increments(): increment of (%d, %d, %d) is not direction %s!!\n\n",
                    ep->increment [0], ep->increment [1], ep->increment [2], clk_dir_name [ep->dir]);
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
                 start [0], start [1], start [2], pos [0], pos [1], pos [2]);

         clk_check_lattice(knot);
         exit(109);
      }
      comp = comp->next;
   }

}



#define NORTH 50
#define EAST  44
#define WEST  41
#define SOUTH 38
#define UP    74
#define DOWN  26

int clk_direction(ivector incr)
{
   // return the direction associated with an increment
   switch ((1 << (incr [0] + 1)) | (1 << (incr [1] + 1 + 2)) | (1 << (incr [2] + 1 + 4)))
   {
      case NORTH:
         return MOVE_NORTH;
         break;
      case EAST:
         return MOVE_EAST;
         break;
      case WEST:
         return MOVE_WEST;
         break;
      case SOUTH:
         return MOVE_SOUTH;
         break;
      case UP:
         return MOVE_UP;
         break;
      case DOWN:
         return MOVE_DOWN;
         break;
   }
   return MOVE_INVALID;
}

int clk_direction(ivector end, ivector start)
{
   ivector incr;
   sub_ivector(incr, end, start);
   return clk_direction(incr);
}

bool clk_check_for_edge_hits(CubicLatticeKnotPtr knot, int dir, ivector test_location)
{
   // we should recentre the knot when an edge hit occurs
   // edge hits don't typically occur except for very small knots
   switch (dir)
   {
      case MOVE_NORTH:
         if (test_location [1] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_SOUTH:
         if (test_location [1] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_EAST:
         if (test_location [0] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_WEST:
         if (test_location [0] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_UP:
         if (test_location [2] > LATTICE_SIZE - 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
      case MOVE_DOWN:
         if (test_location [2] < 2)
         {
            ++knot->edge_hits [dir];
            return true;
         }
         break;
   }

   return false;
}

void clk_get_statistics(CubicLatticeKnotPtr clkp, int &n1, int &n2, int &n3, int &n4)
{
   // gather statistics on local edge configurations with an eye to possibly
   // reordering the cases
   // currently, this function is used only by the KnotPlot driver
   if (!clkp) return;
   n1 = n2 = n3 = n4 = 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         // Case 1: chosen edge is parallel to both adjacent edges,
         // all four moves are +2 moves.

         // directions of adjacent edges same as chosen edge

         if (ep->dir == ep->prev->dir && ep->dir == ep->next->dir)
         {
            ++n1;
         }

            // Case 2: chosen edge is perpendicular to both adjacent edges
            // and those adjacent edges are anti-parallel to each other. 
            // One move is a -2 move, other three are +2 moves.

         else if (perp [ep->dir][ep->next->dir] && anti [ep->prev->dir][ep->next->dir])
         {
            ++n2;
         }

            // Case 3: chosen edge is perpendicular to one adjacent edge
            // and parallel to the other adjacent edge.
            // For this case we know that at least one of adjacent edges is perpendicular to chosen edge
            // because the case of both being parallel has been ruled out above.
            // One move is a 0 move, other three are +2 moves.

            // NOTE: test this before Case 1 above

         else if (ep->dir == ep->prev->dir || ep->dir == ep->next->dir)
         {
            ++n3;
         }

            // Case 4: Only possibility remaining: chosen edge is perpendicular
            // to both adjacent edges, and the adjacent edges are not anti-parallel.
            // Two moves are 0 moves and two moves are +2 moves.

         else
         {
            ++n4;
         }

         ep = ep->next;

      }
      comp = comp->next;
   }
}

bool clkp_always_turns(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (ep->dir == ep->next->dir) return false;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return true;
}

EdgePtr is_tight_clasp(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->dir != ep->prev->dir) return (EdgePtr) NULL;
   if (!perp [ep->dir][ep->next->dir]) return (EdgePtr) NULL;
   if (!anti [ep->prev->prev->dir][ep->next->dir]) return (EdgePtr) NULL;
   ivector loc;
   add_ivector(loc, ep->start, ep->next->increment);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   ivector incr;
   copy_ivector(incr, increment_NEWSUD [kross [ep->dir][ep->next->dir]]);

   ivector loc1;

   add_ivector(loc1, ep->start, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   sub_ivector(loc1, ep->start, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   add_ivector(loc1, loc, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;
   sub_ivector(loc1, loc, incr);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep1 = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(loc, ep1->start))
         {
            if (ep1->dir != ep1->prev->dir) return (EdgePtr) NULL;
            if (!perp [ep1->dir][ep1->next->dir]) return (EdgePtr) NULL;
            if (!anti [ep1->prev->prev->dir][ep1->next->dir]) return (EdgePtr) NULL;
            ivector oloc;
            add_ivector(oloc, ep1->start, ep1->next->increment);
            if (!equal_ivector(oloc, ep->start)) return (EdgePtr) NULL;
            return ep1;
         }
         ep1 = ep1->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL; // should never get here
}

EdgePtr is_parsite_WRONG(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->dir != ep->prev->dir) return (EdgePtr) NULL;

   int xchange_dir;
   ivector loc;
   sub_ivector(loc, ep->start, clkp->loffset);

   switch (ep->dir)
   {
      case MOVE_NORTH:
      case MOVE_SOUTH:
         if (loc [0] > 0)
            xchange_dir = MOVE_EAST;
         else
            xchange_dir = MOVE_WEST;
         break;

      case MOVE_EAST:
      case MOVE_WEST:
         if (loc [1] > 0)
            xchange_dir = MOVE_NORTH;
         else
            xchange_dir = MOVE_SOUTH;
         break;

      default:
         return (EdgePtr) NULL;
   }

   ivector loc2;

   add_ivector(loc, ep->start, increment_NEWSUD [xchange_dir]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [xchange_dir]);
   copy_ivector(loc2, loc);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [ep->dir]);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [opposite [ep->dir]]);
   add_ivector(loc, loc, increment_NEWSUD [opposite [ep->dir]]);
   if (clkp->lattice [lat(loc)] == EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [MOVE_UP]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [opposite [xchange_dir]]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [ep->dir]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [ep->dir]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   add_ivector(loc, loc, increment_NEWSUD [opposite [xchange_dir]]);
   if (clkp->lattice [lat(loc)] != EMPTY) return (EdgePtr) NULL;

   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep1 = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(loc2, ep1->start))
         {
            if (ep1->dir != ep1->prev->dir || ep1->dir != ep->dir) return (EdgePtr) NULL;
            printf("dir found is %s\n", clk_dir_name [xchange_dir]);
            fflush(stdout);
            return ep1;
         }
         ep1 = ep1->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL; // should never get here
}

bool clk_check_path(char *path, CubicLatticeKnotPtr clkp, ivector loc, int value)
{
   ivector loco;
   /*
   printf ("checking path `%s' starting at (%d, %d, %d) for value %d\n",
   path, loc [0] - 128, loc [1] - 128, loc [2] - 128, value); */
   int len = strlen(path);
   copy_ivector(loco, loc);
   for (int i = 0; i < len; i++)
   {
      add_ivector(loco, loco, increment_NEWSUD [(int) (path [i] - '0')]);
      /*
      printf ("moving %d (%s) and looking at (%d, %d, %d)\n", 
         (int) (path [i] - '0'), 
         clk_dir_name [(int) (path [i] - '0')], 
         loco [0]  - 128, loco [1] - 128, loco [2] - 128); */
      if (clkp->lattice [lat(loco)] != value) return false;
   }
   return true;
}

EdgePtr is_parsite(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   if (ep->dir != ep->prev->dir) return (EdgePtr) NULL;

   int xchange_dir;
   ivector loc;
   sub_ivector(loc, ep->start, clkp->loffset);

   switch (ep->dir)
   {
      case MOVE_NORTH:
      case MOVE_SOUTH:
         if (loc [0] > 0)
            xchange_dir = MOVE_EAST;
         else
            xchange_dir = MOVE_WEST;
         break;

      case MOVE_EAST:
      case MOVE_WEST:
         if (loc [1] > 0)
            xchange_dir = MOVE_NORTH;
         else
            xchange_dir = MOVE_SOUTH;
         break;

      default:
         return (EdgePtr) NULL;
   }


   char path [12];
   sprintf(path, "%d%d%d", xchange_dir, ep->dir, ep->dir);
   if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr) NULL;

   sprintf(path, "%d%d%d%d",
           MOVE_UP, ep->dir, xchange_dir, xchange_dir);
   if (!clk_check_path(path, clkp, ep->prev->start, EMPTY)) return (EdgePtr) NULL;

   ivector loc2;
   add_ivector(loc2, ep->start, increment_NEWSUD [xchange_dir]);
   add_ivector(loc2, loc2, increment_NEWSUD [xchange_dir]);
   sprintf(path, "%d%d%d%d",
           MOVE_DOWN, opposite [xchange_dir], opposite [xchange_dir], ep->dir);
   if (!clk_check_path(path, clkp, loc2, EMPTY)) return (EdgePtr) NULL;

   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep1 = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (equal_ivector(loc2, ep1->start))
         {
            if (ep1->dir != ep1->prev->dir || ep1->dir != ep->dir) return (EdgePtr) NULL;
            //printf ("dir found is %s\n", clk_dir_name [xchange_dir]); fflush (stdout);
            return ep1;
         }
         ep1 = ep1->next;
      }
      comp = comp->next;
   }
   return (EdgePtr) NULL; // should never get here
}

/*
void mark_tight_clasp (EdgePtr ep1, EdgePtr ep2) {
  vector3 p1, p2;
  
}
 */
bool has_tight_clasp(CubicLatticeKnotPtr clkp)
{
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (epc = is_tight_clasp(clkp, ep))
            return true;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

void edge_info(char *b, EdgePtr ep)
{
   printf("%s, ID %d, loc %d %d %d, dir %d %d %d : %s, prev %d, next %d, ",
           b, ep->ID,
           ep->start [0] - 128, ep->start [1] - 128, ep->start [2] - 128,
           ep->increment [0], ep->increment [1], ep->increment [2],
           clk_dir_name [ep->dir],
           ep->prev->ID, ep->next->ID);
   if (ep->next->prev == ep)
      printf("Y");
   else
      printf("N");
   if (ep->prev->next == ep)
      printf("Y");
   else
      printf("N");

   ivector a;
   add_ivector(a, ep->start, ep->increment);
   if (equal_ivector(ep->next->start, a))
      printf("Y");
   else
      printf("N");
   add_ivector(a, ep->prev->start, ep->prev->increment);
   if (equal_ivector(ep->start, a))
      printf("Y");
   else
      printf("N");


   printf("\n");
   fflush(stdout);
}

bool has_parsite(CubicLatticeKnotPtr clkp)
{
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (epc = is_parsite(clkp, ep))
            return true;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

void collapse(CubicLatticeKnotPtr clkp, EdgePtr ep)
{
   //  printf ("collapsing edge #%d\n", ep->ID);
   clkp->lattice [lat(ep->prev->start)] = EMPTY;
   clkp->lattice [lat(ep->next->start)] = EMPTY;

   copy_ivector(ep->prev->start, ep->prev->prev->start);
   add_ivector(ep->start, ep->prev->start, ep->prev->increment);

   ComponentCLKPtr comp = (ComponentCLKPtr) ep->comp;
   delete_Edge(clkp, comp, ep->prev->prev);
   delete_Edge(clkp, comp, ep->next);
}

bool bfacf_strand_pass(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   if (clkp->nedges_total < 24) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;

      int shift = rand_integer(0, comp->nedges);
      for (int i = 0; i < shift; i++)
         ep = ep->next;

      for (int i = 0; i < comp->nedges; i++)
      {
         if (epc = is_tight_clasp(clkp, ep))
         {
            collapse(clkp, ep);
            collapse(clkp, epc);
            return true;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

void join_edges(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
   // assumes that start locations are valid
   ep1->next = ep2;
   ep2->prev = ep1;
   sub_ivector(ep1->increment, ep2->start, ep1->start);
   // bbbb
   if (!clk_check_increment(ep1->increment))
   {
      fprintf(stderr, "fails increment check!\n");
      exit(1);
   }
   ep1->dir = clk_direction(ep1->increment);
   clkp->lattice [lat(ep1->start)] = OCCUPIED;
}

void panic_exit_parsite(char *c)
{
   fprintf(stderr, "\n\n *** panick exit!  %s\n", c);
   exit(343);
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp, ComponentCLKPtr comp, EdgePtr ep)
{
   int xchange_dir;
   EdgePtr epc;
   if (!(epc = is_parsite(clkp, ep))) return false;
   if (ep->dir != ep->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 1");
   if (ep->dir != epc->dir) panic_exit_parsite("bfacf_parsite_pass(): 2");
   if (epc->dir != epc->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 3");
   if (ep->dir == MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): 4");
   if (ep->dir == MOVE_DOWN) panic_exit_parsite("bfacf_parsite_pass(): 5");

   edge_info("ep", ep);
   /*
     edge_info ("epp", ep->prev);
    

	

	
 edge_info ("ep", ep);
    edge_info ("epn", ep->next);
    edge_info ("epcp", epc->prev);
    edge_info ("epc", epc);
    edge_info ("epcn", epc->next);
    */

   ivector diff;
   bool found = false;
   sub_ivector(diff, epc->start, ep->start);
   //printf ("diff = %d %d %d\n", diff [0], diff [1], diff [2]);
   for (xchange_dir = 0; xchange_dir < 6; xchange_dir++)
   {
      ivector incr;
      mult_ivector(incr, increment_NEWSUD [xchange_dir], 2);
      //printf ("incr %d = %d %d %d\n", xchange_dir, incr [0], incr [1], incr [2]);
      if (equal_ivector(incr, diff))
      {
         found = true;
         break;
      }
   }
   // bbbb
   if (!found)
   {
      fprintf(stderr, "\n\n *** bfacf_parsite_pass(): xchange_dir not found!\n");
      fflush(stderr);
      exit(934);
   }
   //printf ("dir used is %s\n", clk_dir_name [xchange_dir]); fflush (stdout);

   EdgePtr epn = ep->next;
   EdgePtr epp = ep->prev;
   EdgePtr epcn = epc->next;
   EdgePtr epcp = epc->prev;

   // assign edges to nep[] array as follows:
   //          epp  ep   new edges  epn  epcp  epc  new edges   epcn
   // nep[]:     0   1       2 - 9   10    11   12    13 - 14      15

   // need 10 new edges
   // choose the first 10 unused edges in edge pool
   EdgePtr nep [16];
   int j = 2;
   for (int i = 0; i < 10; i++)
   {
      nep [j] = clkp->edgepool [clkp->nedges_total + i];
      nep [j]->comp = (void *) comp;
      j++;
      if (j == 10) j = 13;
   }
   nep [0] = epp;
   nep [1] = ep;
   nep [10] = epn;
   nep [11] = epcp;
   nep [12] = epc;
   nep [15] = epcn;

   int ep_orig_dir = ep->dir;
   // update positions 
   add_ivector(nep [1]->start, nep [0]->start, increment_NEWSUD [MOVE_UP]);
   add_ivector(nep [2]->start, nep [1]->start, increment_NEWSUD [ep->dir]);
   add_ivector(nep [3]->start, nep [2]->start, increment_NEWSUD [xchange_dir]);
   add_ivector(nep [4]->start, nep [3]->start, increment_NEWSUD [xchange_dir]);
   add_ivector(nep [5]->start, nep [4]->start, increment_NEWSUD [MOVE_DOWN]);
   add_ivector(nep [6]->start, nep [5]->start, increment_NEWSUD [MOVE_DOWN]);
   add_ivector(nep [7]->start, nep [6]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [8]->start, nep [7]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [9]->start, nep [8]->start, increment_NEWSUD [MOVE_UP]);

   add_ivector(nep [12]->start, nep [11]->start, increment_NEWSUD [opposite [xchange_dir]]);
   add_ivector(nep [13]->start, nep [12]->start, increment_NEWSUD [ep->dir]);
   add_ivector(nep [14]->start, nep [13]->start, increment_NEWSUD [ep->dir]);

   // connect edges, adjust adjacency pointers and update lattice
   for (int i = 0; i < 15; i++)
   {
      if (i == 10) continue;
      join_edges(clkp, nep [i], nep [i + 1]);
   }

   // increase by 10 number of edges in comp and in knot
   comp->nedges += 10;
   clkp->nedges_total += 10;

   return true;
}

bool bfacf_parsite_pass(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   if (clkp->nedges_total < 24) return false;
   if (clkp->ncomps > 1) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;

      int shift = rand_integer(0, comp->nedges);
      //printf ("shifting by %d\n", shift);
      for (int i = 0; i < shift; i++)
         ep = ep->next;

      for (int nn = 0; nn < comp->nedges; nn++)
      {
         if (bfacf_parsite_pass(clkp, comp, ep)) return true;
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

bool bfacf_parsite_pass_WRONG(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   if (clkp->nedges_total < 24) return false;
   if (clkp->ncomps > 1) return false;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      int xchange_dir;
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;

      int shift = rand_integer(0, comp->nedges);
      for (int i = 0; i < shift; i++)
         ep = ep->next;

      for (int nn = 0; nn < comp->nedges; nn++)
      {
         if (epc = is_parsite(clkp, ep))
         {
            if (ep->dir != ep->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 1");
            if (ep->dir != epc->dir) panic_exit_parsite("bfacf_parsite_pass(): 2");
            if (epc->dir != epc->prev->dir) panic_exit_parsite("bfacf_parsite_pass(): 3");
            if (ep->dir == MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): 4");
            if (ep->dir == MOVE_DOWN) panic_exit_parsite("bfacf_parsite_pass(): 5");

            edge_info("epp", ep->prev);
            edge_info("ep", ep);
            edge_info("epn", ep->next);
            edge_info("epcp", epc->prev);
            edge_info("epc", epc);
            edge_info("epcn", epc->next);

            ivector diff;
            bool found = false;
            sub_ivector(diff, epc->start, ep->start);
            //printf ("diff = %d %d %d\n", diff [0], diff [1], diff [2]);
            for (xchange_dir = 0; xchange_dir < 6; xchange_dir++)
            {
               ivector incr;
               mult_ivector(incr, increment_NEWSUD [xchange_dir], 2);
               //printf ("incr %d = %d %d %d\n", xchange_dir, incr [0], incr [1], incr [2]);
               if (equal_ivector(incr, diff))
               {
                  found = true;
                  break;
               }
            }
            // bbbb
            if (!found)
            {
               fprintf(stderr, "\n\n *** bfacf_parsite_pass(): xchange_dir not found!\n");
               fflush(stderr);
               exit(934);
            }
            printf("dir used is %s\n", clk_dir_name [xchange_dir]);
            fflush(stdout);

            EdgePtr epn = ep->next;
            //EdgePtr epnn = ep->next->next;
            EdgePtr epcp = epc->prev;

            // need 6 new edges
            // choose the first 6 unused edges in edge pool
            EdgePtr nep [6];
            for (int i = 0; i < 6; i++)
            {
               nep [i] = clkp->edgepool [clkp->nedges_total + i];
               nep [i]->comp = (void *) comp;
            }

            int ep_orig_dir = ep->dir;
            // update positions 
            add_ivector(nep [0]->start, ep->start, increment_NEWSUD [xchange_dir]);
            add_ivector(nep [1]->start, epcp->start, increment_NEWSUD [MOVE_UP]);
            add_ivector(nep [2]->start, nep [1]->start, increment_NEWSUD [opposite [xchange_dir]]);
            add_ivector(nep [3]->start, nep [2]->start, increment_NEWSUD [ep->dir]);
            add_ivector(nep [4]->start, nep [3]->start, increment_NEWSUD [ep->dir]);
            add_ivector(nep [5]->start, nep [4]->start, increment_NEWSUD [opposite [xchange_dir]]);

            // connect edges, adjust adjacency pointers and update lattice
            join_edges(clkp, ep, nep [0]);
            join_edges(clkp, nep [0], epc);
            join_edges(clkp, epcp, nep [1]);
            join_edges(clkp, nep [1], nep [2]);
            join_edges(clkp, nep [2], nep [3]);
            join_edges(clkp, nep [3], nep [4]);
            join_edges(clkp, nep [4], nep [5]);
            join_edges(clkp, nep [5], epn);

            if (ep->dir != xchange_dir) panic_exit_parsite("bfacf_parsite_pass(): A");
            if (nep [0]->dir != xchange_dir) panic_exit_parsite("bfacf_parsite_pass(): B");
            if (epcp->dir != MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): C");
            if (nep [1]->dir != opposite [xchange_dir]) panic_exit_parsite("bfacf_parsite_pass(): D");
            if (nep [2]->dir != ep_orig_dir) panic_exit_parsite("bfacf_parsite_pass(): E");
            if (nep [3]->dir != ep_orig_dir) panic_exit_parsite("bfacf_parsite_pass(): F");
            if (nep [4]->dir != opposite [xchange_dir]) panic_exit_parsite("bfacf_parsite_pass(): G");
            if (nep [5]->dir != MOVE_DOWN) panic_exit_parsite("bfacf_parsite_pass(): H");
            if (xchange_dir == MOVE_UP) panic_exit_parsite("bfacf_parsite_pass(): I");
            // increase by 6 number of edges in comp and in knot
            comp->nedges += 6;
            clkp->nedges_total += 6;

            edge_info("ep", ep);
            edge_info("epc", epc);
            edge_info("epcp", epcp);
            edge_info("epn", epn);
            for (int i = 0; i < 6; i++)
            {
               char name [8];
               sprintf(name, "nep %d", i);
               edge_info(name, nep [i]);
            }

            clk_check_increments(clkp);
            return true;
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return false;
}

int numb_tight_clasps(CubicLatticeKnotPtr clkp)
{
   return numb_tight_clasps(clkp, NULL);
}

int numb_tight_clasps(CubicLatticeKnotPtr clkp, void mark(ivector, ivector))
{
   int kount = 0;
   if (!clkp) return 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (epc = is_tight_clasp(clkp, ep))
         {
            if (ep->ID > epc->ID)
            {
               ++kount;
               if (mark)
                  mark(ep->start, epc->start);
            }
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return kount;
}

int numb_parsites(CubicLatticeKnotPtr clkp)
{
   return numb_parsites(clkp, NULL);
}

int numb_parsites(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr))
{
   int kount = 0;
   if (!clkp) return 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      EdgePtr epc;
      for (int i = 0; i < comp->nedges; i++)
      {
         if (epc = is_parsite(clkp, ep))
         {
            ++kount;
            if (mark)
               mark(ep, epc);
         }
         ep = ep->next;
      }
      comp = comp->next;
   }
   return kount;
}

void init_lattice(CubicLatticeKnotPtr clkp, ivector mid)
{
   // offset the knot so that the point `mid' lies at at 
   // (LATTICE_SIZE/2, LATTICE_SIZE/2, LATTICE_SIZE/2)
   set_ivector(clkp->loffset, LATTICE_SIZE / 2, LATTICE_SIZE / 2, LATTICE_SIZE / 2);
   sub_ivector(clkp->loffset, clkp->loffset, mid);
   set_ivector(clkp->max_range, LATTICE_SIZE - 2, LATTICE_SIZE - 2, LATTICE_SIZE - 2);
   ComponentCLKPtr comp = clkp->fcomp;

   while (comp)
   {
      // now go thru all the edges, adding the loffset and marking as OCCUPIED
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         add_ivector(ep->start, ep->start, clkp->loffset); // we will subtract this offset when saving the knot
         clkp->lattice [lat(ep->start)] = OCCUPIED;
         ep = ep->next;
      }
      comp = comp->next;
   }
}

//#define RGS_PARANOIA

bool recentre_knot_in_lattice(CubicLatticeKnotPtr knot)
{ // returns true if knot can't be recentered
   if (!knot) return true;

   ComponentCLKPtr comp;
   ivector min, max, mid;
   ivector diff;
   EdgePtr ep;

   clk_get_extent(min, max, knot);
   if (max [0] - min [0] > knot->max_range [0]) return true;
   if (max [1] - min [1] > knot->max_range [1]) return true;
   if (max [2] - min [2] > knot->max_range [2]) return true;

   // first clear the old path
   clear_lattice(knot);
   ep = knot->fcomp->first_edge;

#ifdef RGS_PARANOIA
   sub_ivector(diff, ep->start, knot->loffset);
   printf("before: %d %d %d   %d %d %d   %d %d %d\n",
           ep->start [0], ep->start [1], ep->start [2],
           knot->loffset [0], knot->loffset [1], knot->loffset [2],
           diff [0], diff [1], diff [2]);
   fflush(stdout);
#endif

   comp = knot->fcomp;
   while (comp)
   {
      ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         // recover the start locations for the edges in actual 3D space
         sub_ivector(ep->start, ep->start, knot->loffset);
         ep = ep->next;
      }
      comp = comp->next;
   }

   midpoint(mid, max, min);

   init_lattice(knot, mid);

#ifdef RGS_PARANOIA
   ep = knot->fcomp->first_edge;
   sub_ivector(diff, ep->start, knot->loffset);
   printf("after: %d %d %d   %d %d %d   %d %d %d\n\n",
           ep->start [0], ep->start [1], ep->start [2],
           knot->loffset [0], knot->loffset [1], knot->loffset [2],
           diff [0], diff [1], diff [2]);
   fflush(stdout);
#endif

   return false;
}

// currently, set_lattice(), fill_lattice(), invert_lattice(), dilate_lattice(), and clear_lattice()
// are used only by the KnotPlot driver to implement topological obstructions and cavities

void set_lattice(CubicLatticeKnotPtr clkp, char value, ivector min, ivector max)
{
   if (!clkp) return;
   //clkp->auto_recentre = false;
   ivector offset_min, offset_max;
   add_ivector(offset_min, min, clkp->loffset);
   add_ivector(offset_max, max, clkp->loffset);
   for (int i = 0; i < 3; i++)
   {
      CLAMP(offset_min [i], 0, LATTICE_SIZE - 1);
      CLAMP(offset_max [i], 0, LATTICE_SIZE - 1);
   }
   ivector a;

   for (a [0] = offset_min [0]; a [0] <= offset_max [0]; a [0]++)
      for (a [1] = offset_min [1]; a [1] <= offset_max [1]; a [1]++)
         for (a [2] = offset_min [2]; a [2] <= offset_max [2]; a [2]++)
            clkp->lattice [lat(a)] = value;
}

void fill_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max)
{
   set_lattice(clkp, (char) OCCUPIED, min, max);
}

void invert_lattice(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return;
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (clkp->lattice [i] == OCCUPIED)
         clkp->lattice [i] = EMPTY;
      else
         clkp->lattice [i] = OCCUPIED;
   }
}

void dilate_lattice(CubicLatticeKnotPtr clkp)
{
   // expand occupied region to occupy neighboring cells
   if (!clkp) return;
   ivector pos;

   for (pos [2] = 1; pos [2] < LATTICE_SIZE - 1; pos [2]++)
   { // unnecessesary brace to make VC++ indenter happy
      for (pos [1] = 1; pos [1] < LATTICE_SIZE - 1; pos [1]++)
      {
         for (pos [0] = 1; pos [0] < LATTICE_SIZE - 1; pos [0]++)
         {
            if (clkp->lattice [lat(pos)] == OCCUPIED)
            {
               for (int i = 0; i < 6; i++)
               {
                  ivector pos2;
                  add_ivector(pos2, pos, increment_NEWSUD [i]);
                  if (clkp->lattice [lat(pos2)] == EMPTY)
                  {
                     clkp->lattice [lat(pos2)] = NEXTOCCUPIED;
                  }
               }
            }
         }
      }
   }
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
   {
      if (clkp->lattice [i] == NEXTOCCUPIED)
         clkp->lattice [i] = OCCUPIED;
   }
}

int bfacf_info_filled(CubicLatticeKnotPtr knot)
{
   int filled = 0;
   for (int i = 0; i < LATTICE_TOTAL_SIZE; i++)
      if (knot->lattice [i] == OCCUPIED)
         ++filled;
   return filled;
}

void set_lattice(CubicLatticeKnotPtr clkp, char value)
{
   if (!clkp) return;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      ivector a;
      for (int i = 0; i < comp->nedges; i++)
      {
         clkp->lattice [lat(ep->start)] = value;
         ep = ep->next;
      }
      if (comp->flags & COMPONENT_CLK_FLAG_OPEN)
      {
         add_ivector(a, ep->start, ep->increment);
         clkp->lattice [lat(a)] = value;
      }
      comp = comp->next;
   }
}

void fill_lattice(CubicLatticeKnotPtr clkp)
{
   set_lattice(clkp, OCCUPIED);
}

void clear_lattice(CubicLatticeKnotPtr clkp)
{
   set_lattice(clkp, EMPTY);
}

void clear_lattice(CubicLatticeKnotPtr clkp, ivector min, ivector max)
{
   set_lattice(clkp, (char) EMPTY, min, max);
}

void bfacf_set_lattice_sphere(CubicLatticeKnotPtr clkp, int what, double radius, double xcen, double ycen, double zcen, char *blurt)
{
#ifdef INCLUDED_FROM_KnotPlot 
   if (!clkp) return;
   char action [14];
   if (blurt)
   {
      if (what == OCCUPIED)
         sprintf(action, "filling");
      else
         sprintf(action, "clearing");
      sprintf(blurt, "%s lattice with sphere of radius %.2f centered at (%.2f, %.2f, %.2f)\n", action, radius, xcen, ycen, zcen);
   }

   ivector pos;
   vector3 cen;
   set_vector(cen, xcen + (double) clkp->loffset [0], ycen + (double) clkp->loffset [1], zcen + (double) clkp->loffset [2]);
   double r2 = radius * radius;

   for (pos [2] = 1; pos [2] < LATTICE_SIZE - 1; pos [2]++)
   { // unnecessesary brace to make VC++ indenter happy
      for (pos [1] = 1; pos [1] < LATTICE_SIZE - 1; pos [1]++)
      {
         for (pos [0] = 1; pos [0] < LATTICE_SIZE - 1; pos [0]++)
         {
            vector3 rpos;
            copy_vector(rpos, pos);
            if (diff_distance_SQ(rpos, cen) < r2)
               clkp->lattice [lat(pos)] = what;
         }
      }
   }
#endif
}

/*
void compute_energy (CubicLatticeKnotPtr clkp) {
  double energy = 0.0;
  EdgePtr ep1 = clkp->first_edge;
  
}
 */

void get_NEWS(char *news, CubicLatticeKnotPtr clkp)
{
   if (!clkp) return;
   if (clkp->ncomps > 1)
   {
      strcpy(news, "unknown");
      return;
   }
   ComponentCLKPtr comp = clkp->fcomp;
   EdgePtr ep = comp->first_edge;

   for (int i = 0; i < comp->nedges; i++)
   {
      news [i] = clk_dir_name [ep->dir][0];
      ep = ep->next;
   }
   news [comp->nedges] = (char) 0;
}

void copy_CLKP_to_array(ivector *coord, int &nedges, CubicLatticeKnotPtr clkp)
{
   if (!clkp) return;
   ComponentCLKPtr comp = clkp->fcomp;
   nedges = 0;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      nedges += comp->nedges;
      for (int i = 0; i < comp->nedges; i++)
      {
         sub_ivector(coord [i], ep->start, clkp->loffset);
         ep = ep->next;
      }
      comp = comp->next;
   }
}

void print_zknown(bool verbose)
{
   printf("z-values are known for the following lengths:");
   for (int i = 0; i < NUM_KNOWN_LENGTHS; i++)
      printf(" %d", known_length [i]);
   printf("\nNote that for a given length, the z-value might not be known for some knot types.\n");
}

bool findz(int target_length, char *knot2, double &z)
{
   char knot [8];
   int len = (int) strlen(knot2);
   if (len > 5) return false;
   double *zvalue;
   strcpy(knot, knot2);
   if (knot [len - 1] == 's') knot [len - 1] = (char) 0;
   switch (target_length)
   {
      case 56:
         zvalue = zvalue_56;
         break;
      case 60:
         zvalue = zvalue_60;
         break;
      case 65:
         zvalue = zvalue_65;
         break;
      case 70:
         zvalue = zvalue_70;
         break;
      case 75:
         zvalue = zvalue_75;
         break;
      case 80:
         zvalue = zvalue_80;
         break;
      case 85:
         zvalue = zvalue_85;
         break;
      case 90:
         zvalue = zvalue_90;
         break;
      case 95:
         zvalue = zvalue_95;
         break;
      case 100:
         zvalue = zvalue_100;
         break;
      case 105:
         zvalue = zvalue_105;
         break;
      case 110:
         zvalue = zvalue_110;
         break;
      case 115:
         zvalue = zvalue_115;
         break;
      case 120:
         zvalue = zvalue_120;
         break;
      case 124:
         zvalue = zvalue_124;
         break;
      case 125:
         zvalue = zvalue_125;
         break;
      case 130:
         zvalue = zvalue_130;
         break;
      case 135:
         zvalue = zvalue_135;
         break;
      case 140:
         zvalue = zvalue_140;
         break;
      case 145:
         zvalue = zvalue_145;
         break;
      case 150:
         zvalue = zvalue_150;
         break;
      case 155:
         zvalue = zvalue_155;
         break;
      case 160:
         zvalue = zvalue_160;
         break;
      case 165:
         zvalue = zvalue_165;
         break;
      case 170:
         zvalue = zvalue_170;
         break;
      case 175:
         zvalue = zvalue_175;
         break;
      case 180:
         zvalue = zvalue_180;
         break;
      case 185:
         zvalue = zvalue_185;
         break;
      case 190:
         zvalue = zvalue_190;
         break;
      case 195:
         zvalue = zvalue_195;
         break;
      case 200:
         zvalue = zvalue_200;
         break;
      case 300:
         zvalue = zvalue_300;
         break;
      default:
         return false;
   }
   for (int i = 0; i < NUM_KNOTS; i++)
   {
      if (!strcmp(knotname [i], knot))
      {
         z = zvalue [i];
         if (z > 1.0) return false;
         return true;
      }
   }

   return false;
}

char *get_known_lengths(char *knot)
{
   char *s = (char *) calloc(NUM_KNOWN_LENGTHS * 4 + 2, sizeof (char));
   char buff [8];
   double dummy;
   for (int l = 0; l < NUM_KNOWN_LENGTHS; l++)
   {
      if (findz(known_length [l], knot, dummy))
      {
         sprintf(buff, "%d ", known_length [l]);
         strcat(s, buff);
      }
   }
   return s;
}

static double kappa = DEFAULT_KAPPA; // a value used in the paper
static double nu = -0.26; // a value used in the paper
static double A = DEFAULT_A; // change 0.01 to 0.1 and 1.0   

#define ABS(X)      ((X) < 0 ? -(X) : (X))


double T = DEFAULT_KELVINS; // approximately body temperature

void energy_set_nu(double value)
{
   nu = value;
}

void energy_set_T(double value)
{
   T = value;
}

void energy_set_kappa(double value)
{
   kappa = value;
}

void energy_set_A(double value)
{
   A = value;
}

void energy_blurt_values(char *s)
{
   sprintf(s, "A = %f, kappa = %f, nu = %f", A, kappa, nu);
}

double energy(ivector *currentKnot, int length)
{
   int i, j;
   int distance; // distance is back to being an int 
   double total_energy = 0.0;

   int last = length - 1;

   for (i = 0; i < length - 1; i++)
   {
      for (j = i + 2; j < last; j++)
      { // RGS: i + 2 because we don't want adjacent vertices
         distance = ABS(currentKnot [i][0] - currentKnot [j][0]) +
                 ABS(currentKnot [i][1] - currentKnot [j][1]) +
                 ABS(currentKnot [i][2] - currentKnot [j][2]);

         total_energy += A * exp(-kappa * (double) distance) / distance;

         if (distance == 1)
            total_energy += nu;

      }
      last = length;
   }

   return total_energy;
}

static double *expd;

void bfacf_initialize_expd()
{
   printf("initializing with A = %f and kappa = %f\n", A, kappa);
   expd = (double *) calloc(3 * LATTICE_SIZE, sizeof (double));
   for (int distance = 1; distance < 3 * LATTICE_SIZE; distance++)
      expd [distance] = A * exp(-kappa * (double) distance) / (double) distance;
}

double energy(ivector location, CubicLatticeKnotPtr clkp, EdgePtr beg, EdgePtr end)
{
   // compare `location' to the starting location of each edge from `beg' to `end'

   EdgePtr ep = beg;
   double energy = 0.0;

   while (true)
   {
      int distance =
              ABS(location [0] - ep->start [0]) +
              ABS(location [1] - ep->start [1]) +
              ABS(location [2] - ep->start [2]);

      energy += expd [distance];

      if (distance == 1)
         energy += nu;

      if (ep == end) break;
      ep = ep->next;
   }

   return energy;
}

double energy(ivector locationA, ivector locationB)
{
   double energy = 0.0;

   int distance =
           ABS(locationA [0] - locationB [0]) +
           ABS(locationA [1] - locationB [1]) +
           ABS(locationA [2] - locationB [2]);

   energy = expd [distance];

   if (distance == 1)
      energy += nu;

   return energy;
}

void energy_report(CubicLatticeKnotPtr clkp, char *report)
{
   if (!clkp) return;
   ComponentCLKPtr comp = clkp->fcomp;
   ivector *loc = (ivector *) calloc(comp->nedges, sizeof (ivector));
   EdgePtr ep = comp->first_edge;
   for (int i = 0; i < comp->nedges; i++)
   {
      copy_ivector(loc [i], ep->start);
      ep = ep->next;
   }
   sprintf(report, "A = %f, kappa = %f, nu = %f, E = %f\n", A, kappa, nu,
           energy(loc, comp->nedges));
   free(loc);
}

bool MMC_swap(double &prob, double zi, int lengthi, double zip1, int lengthip1)
{
   if (lengthip1 < lengthi) return true;
   //return false;
   // chain i has z-value `zi' and length `lengthi'
   // chain i+1 has z-value `zip1' and length `lengthip1'
   //prob = pow (zi, (double) (lengthip1 - lengthi)) * pow (zip1, (double) (lengthi - lengthip1));
   prob = pow((zi / zip1), (double) (lengthip1 - lengthi)); // should precomute this into a table, or into tables
   if (rand_integer(0, 100000) == 102)
   {
      printf("%.4f ", prob);
      fflush(stdout);
   }

   if (rand_uniform() < prob)
      return true;
   else
      return false;
}

void MMC_prob_report(char *knot)
{
   printf("knot %s\n\n", knot);
   for (int l1 = 0; l1 < NUM_KNOWN_LENGTHS; l1++)
   {
      double z1, z2;
      if (!findz(known_length [l1], knot, z1)) continue;
      printf("length %d at z = %f\n", known_length [l1], z1);
      for (int l2 = l1 + 1; l2 < NUM_KNOWN_LENGTHS; l2++)
      {
         if (!findz(known_length [l2], knot, z2)) continue;
         double p;
         MMC_swap(p, z1, l1, z2, l2);
         printf("%d %f %.2f\n", known_length [l2], z2, 100.0 * p);
      }
      printf("\n");
   }
}

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

//qProb *q_prob = (qProb *) NULL;
//int q_prob_size = 0;

void init_probabilities_q(CubicLatticeKnotPtr clkp, int S, double z, double q, bool blurt)
{
   /*
   static int last_S = 0;
   static double last_z = 102.2;
   static double last_q = 102.2;
   if (last_S == S && last_z == z && last_q == q && q_prob) return;
   if (q_prob && S > last_S) {
     free (q_prob);
     q_prob = (qProb *) NULL;
     if (blurt) 
       printf ("freeing previous q_prob\n");
   }

   q_prob_size = last_S = S;
   last_z = z;
   last_q = q;
    */

   if (blurt)
   {
      printf("*** initializing q_prob with S = %d, z = %f, q = %f\n", S, z, q);
      fflush(stdout);
   }

   clkp->q_prob_size = S;
   if (!clkp->q_prob)
      clkp->q_prob = (qProb *) calloc(clkp->q_prob_size, sizeof (qProb));
   for (int ni = 4; ni < clkp->q_prob_size; ni++)
   {
      double n = (double) ni;
      clkp->q_prob [ni].p_plus2 = pow(n + 2, q - 1) * (z * z) / (pow(n, q - 1) + 3.0 * pow(n + 2, q - 1) * z * z);
      //q_prob [ni].p_0      = (1.0 + z * z) / (2.0 * (1.0 + 3.0 * z * z));
      clkp->q_prob [ni].p_minus2 = pow(n - 2, q - 1) / (pow(n - 2, q - 1) + 3.0 * pow(n, q - 1) * z * z);
      clkp->q_prob [ni].p_0 = (clkp->q_prob [ni].p_plus2 + clkp->q_prob [ni].p_minus2) / 2.0;

      clkp->q_prob [ni].p_4p2 = 4.0 * clkp->q_prob [ni].p_plus2;
      clkp->q_prob [ni].p_03p2 = clkp->q_prob [ni].p_0 + 3.0 * clkp->q_prob [ni].p_plus2;
      clkp->q_prob [ni].p_m23p2 = clkp->q_prob [ni].p_minus2 + 3.0 * clkp->q_prob [ni].p_plus2;
      clkp->q_prob [ni].p_2p0 = 2.0 * clkp->q_prob [ni].p_0;
      clkp->q_prob [ni].p_2p02p2 = 2.0 * clkp->q_prob [ni].p_0 + 2.0 * clkp->q_prob [ni].p_plus2;
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

//  
#ifdef INCLUDED_FROM_KnotPlot

void showem_obsolete(CubicLatticeKnotPtr clkp, int kount)
{
   extern void decode_line(char *);
   char com [44];
   decode_line("txt delete group 999");
   decode_line("txt group 999");
   decode_line("txt font 4");
   decode_line("txt fill off");
   decode_line("mode l");
   decode_line("show none");
   extern void output_it(char *);
   decode_line("scen dscal .2");
   sprintf(com, "nfrozen is %d\n", clkp->nfrozen);
   output_it(com);
   for (int i = 0; i < clkp->nfrozen + kount; i++)
   {
      vector3 v;
      EdgePtr ep = clkp->edgepool [i];
      ivector a, b;
      sub_ivector(a, ep->start, clkp->loffset);
      sub_ivector(b, ep->next->start, clkp->loffset);
      midpoint(v, a, b);
      sprintf(com, "txt %f %f %f e%d", v [0], v [1], v[2], ep->locpool);
      decode_line(com);
   }
   decode_line("txt group 0");
}
#endif

void pexit(char *s)
{
   fprintf(stderr, "%s\n", s);
   fflush(stderr);
   exit(102);
}

bool clk_recombo_direct = true;

void clk_fix_incr(EdgePtr ep)
{
   sub_ivector(ep->increment, ep->next->start, ep->start);
   ep->dir = clk_direction(ep->increment);
   if (ep->dir == MOVE_INVALID)
      pexit("invalid increment detected");
}

bool perform_recombination_inverted(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
   if (ep1->comp != ep2->comp) return false;
   ivector test;
   sub_ivector(test, ep1->start, ep2->start);

   // the following section can be removed once we use pairs
   if (!clk_check_increment(test))
      pexit("impossible situation perform_recombination_inverted () A");
   if (ep1->next == ep2 || ep1->prev == ep2)
      pexit("impossible situation perform_recombination_inverted () B");


   EdgePtr ep1n = ep1->next;
   EdgePtr ep2n = ep2->next;
   EdgePtr ep1nn = ep1n->next;
   EdgePtr ep2p = ep2->prev;

   EdgePtr epp = ep1n;
   EdgePtr ep = ep1n->next;

   int kount = 0;

   while (ep != ep2)
   {
      EdgePtr epn = ep->next;
      ep->next = epp;
      ep->prev = epn;
      epp = ep;
      ep = epn;
      if (++kount > 16008)
         pexit("infinite loop A");
   }

   ep1->next = ep2;
   ep2->prev = ep1;
   ep2->next = ep2p;
   ep1n->prev = ep1nn;
   ep1n->next = ep2n;
   ep2n->prev = ep1n;

   /*
   ivector ep1start;
   ivector ep2start;
   sub_ivector (ep1start, ep1->start, clkp->loffset);
   sub_ivector (ep2start, ep2->start, clkp->loffset);
   printf ("%d %d %d, %d %d %d, kount = %d, nedges = %d\n", 
      ep1start [0], ep1start [1], ep1start [2], 
      ep2start [0], ep2start [1], ep2start [2], 
      kount, clkp->fcomp->nedges);*/

   kount = 0;
   ep = ep1;
   do
   {
      clk_fix_incr(ep);
      ep = ep->next;
      if (++kount > 16008)
         pexit("infinite loop B");
   }
   while (ep != ep2n);

   clk_validate(clkp, "end of perform_recombination_inverted()");
   return true;
}

static int null_bug_a, null_bug_b, null_bug_c;

bool perform_recombination(CubicLatticeKnotPtr clkp, EdgePtr ep1, EdgePtr ep2)
{
   if (!ep1 || !ep2)
   {
#ifdef INCLUDED_FROM_KnotPlot
      complain("!perform_recombination (): null pointer");
#endif
      fprintf(stderr, " *** perform_recombination (): null pointer");
      if (!ep1) fprintf(stderr, "ep1 is null ");
      if (!ep2) fprintf(stderr, "ep2 is null ");
      fflush(stderr);

      return false;
   }
   // following is needed to prevent edges like 29 and 31 in 9jun10a.k from being recomboed
   if (ep1->next == ep2->prev || ep1->prev == ep2->next) return false;

   if (ep1->dir == ep2->dir)
      return perform_recombination_inverted(clkp, ep1, ep2);

   // do some paranoia checking (for now)
   if (!anti [ep1->dir][ep2->dir])
      pexit("edges are not anti-parallel as expected!");


   ivector incr;
   sub_ivector(incr, ep1->next->start, ep2->start);
   if (!clk_check_increment(incr))
      pexit("invalid increment found!");

   ComponentCLKPtr comp1 = (ComponentCLKPtr) ep1->comp;
   ComponentCLKPtr comp2 = (ComponentCLKPtr) ep2->comp;

   if (comp1 == comp2 && comp1->nedges < 5) return false;
   if (clk_recombo_limit && comp1 == comp2 && clkp->ncomps > 1) return false;

   EdgePtr ep1n = ep1->next;
   EdgePtr ep1p = ep1->prev;
   EdgePtr ep2n = ep2->next;
   EdgePtr ep2p = ep2->prev;

   // update directions of edges
   sub_ivector(ep1->increment, ep2->next->start, ep1->start);
   ep1->dir = clk_direction(ep1->increment);
   sub_ivector(ep2->increment, ep1->next->start, ep2->start);
   ep2->dir = clk_direction(ep2->increment);

   // update how things are connected
   ep1->next = ep2n;
   ep2->next = ep1n;
   ep1n->prev = ep2;
   ep2n->prev = ep1;

   if (comp1 == comp2)
   {
      comp2 = (ComponentCLKPtr) calloc(1, sizeof (ComponentCLK));
      comp2->minedges = clk_minedges;
      comp2->maxedges = clk_maxedges;
      memcpy(comp2, comp1, sizeof (ComponentCLK));
      comp1->first_edge = ep1;
      comp1->last_edge = ep1->prev;
      comp2->first_edge = ep2;
      comp2->last_edge = ep2->prev;
      comp1->nedges = clk_count_edges(ep1);
      comp2->nedges = clk_count_edges(ep2);
#ifdef MACOSX
      if (comp1->nedges < 4 || comp2->nedges < 4) system("say oh dear, not again &");
#endif
      // add new component comp2 to end of linked list
      comp2->prev = clkp->lcomp;
      comp2->next = (ComponentCLKPtr) NULL;
      clkp->lcomp->next = comp2;
      clkp->lcomp = comp2;
      clkp->ncomps++;
   }
   else
   {
      comp1->first_edge = ep1;
      comp1->last_edge = ep1->prev;
      comp1->nedges = clk_count_edges(ep1);

      // remove comp2 from linked list
      if (comp2->prev)
         comp2->prev->next = comp2->next;
      else
         clkp->fcomp = comp2->next;

      if (comp2->next)
         comp2->next->prev = comp2->prev;
      else
         clkp->lcomp = comp2->prev;

      free(comp2);

      clkp->ncomps--;
   }

   // recompute IDs
   ComponentCLKPtr comp = clkp->fcomp;
   int ID = 0;
   while (comp)
   {
      comp->ID = ID++;
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         ep->comp = (void *) comp;
         ep = ep->next;
      }
      comp = comp->next;
   }

   //  validate (clkp, "end of perform_recombination()");

   return true;
}

//static int recombo_pair_index = 0;
//static EdgePairPtr recombo_pair = (EdgePairPtr) NULL;
//#define DEFAULT_RECOMBO_MAX 108
//static int recombo_max = 0;
//static bool recomboInitialized = false;
//
//int getNumberAvailableRecomboPairs(void)
//{
//   return recombo_pair_index;
//}
//
//void initRecombo(void)
//{
//   recombo_pair_index = 0;
//   if (!recomboInitialized)
//   {
//      recombo_pair = (EdgePairPtr) calloc(DEFAULT_RECOMBO_MAX, sizeof (EdgePair));
//      recombo_max = DEFAULT_RECOMBO_MAX;
//      recomboInitialized = true;
//   }
//}
//
//EdgePairPtr getAvailableRecomboPairs(void)
//{
//   return recombo_pair;
//}
//
//void recombo_record_pair(EdgePtr ep1, EdgePtr ep2)
//{
//   if (!ep1)
//   {
//#ifdef INCLUDED_FROM_KnotPlot
//      complain("!recombo_record_pair(): ep1 NULL POINTER!");
//#endif
//      fprintf(stderr, " *** recombo_record_pair(): ep1 NULL POINTER! recombo_max = %d, recombo_pair_index = %d", recombo_max, recombo_pair_index);
//      return;
//   }
//   if (!ep2)
//   {
//#ifdef INCLUDED_FROM_KnotPlot
//      complain("!recombo_record_pair(): ep2 NULL POINTER!");
//#endif
//      fprintf(stderr, " *** recombo_record_pair(): ep2 NULL POINTER! recombo_max = %d, recombo_pair_index = %d", recombo_max, recombo_pair_index);
//      return;
//   }
//   if (recombo_pair_index == recombo_max)
//   {
//      int new_size = (140 * recombo_max) / 100 + 108;
//      EdgePairPtr epp = (EdgePairPtr) calloc(new_size, sizeof (EdgePair));
//      if (recombo_pair)
//      {
//         fprintf(stderr, "recombo_record_pair(): %d\n", new_size);
//         for (int i = 0; i < recombo_max; i++)
//         {
//            epp [i].ep1 = recombo_pair [i].ep1;
//            epp [i].ep2 = recombo_pair [i].ep2;
//         }
//         free(recombo_pair);
//      }
//      recombo_pair = epp;
//      recombo_max = new_size;
//      ++recombo_pair_index;
//      return;
//   }
//   if (ep1 == ep2)
//      pexit("recombo_record_pair()");
//   recombo_pair [recombo_pair_index].ep1 = ep1;
//   recombo_pair [recombo_pair_index].ep2 = ep2;
//   ++recombo_pair_index;
//}
//
//int bfacf_perform_recombination(CubicLatticeKnotPtr clkp)
//{
//   recombo_pair_index = 0;
//   // find all the pairs of edges involved in recombo sites and record them
//   int numb = bfacf_perform_recombination(clkp, recombo_record_pair);
//   if (numb < 1) return 0; // no recombo sites
//
//   // choose a pair at random
//   int pair = rand_integer(0, numb);
//   if (pair < 0 || pair >= numb)
//   {
//#ifdef INCLUDED_FROM_KnotPlot
//      complain("!bad selection!");
//#endif
//      fprintf(stderr, "pair = 0, numb = %d\n", pair, numb);
//      return 0;
//   }
//   perform_recombination(clkp, recombo_pair [pair].ep1, recombo_pair [pair].ep2);
//
//   return numb;
//}

void clk_fill_path_alt(CubicLatticeKnotPtr clkp, bool fill)
{
   // fill knot path with edge pool location 
   int mask = 0xffffff;
   if (!fill)
      mask = 0;
   ComponentCLKPtr comp = clkp->fcomp;
   while (comp)
   {
      EdgePtr ep = comp->first_edge;
      for (int i = 0; i < comp->nedges; i++)
      {
         clkp->alt_lattice [lat(ep->start)] = (ep->locpool + 1) & mask;
         ep = ep->next;
      }
      comp = comp->next;
   }
}

bool clk_allocation_alt_lattice(CubicLatticeKnotPtr clkp)
{
   if (!clkp) return false;
   if (!clkp->alt_lattice)
      clkp->alt_lattice = (int *) calloc(LATTICE_TOTAL_SIZE, sizeof (int));
   return true;
}

int bfacf_perform_recombination(CubicLatticeKnotPtr clkp, void mark(EdgePtr, EdgePtr))
{
   if (!clk_allocation_alt_lattice(clkp)) return 0;


   //  validate (clkp, "start of A");

   // fill knot path with edge pool location 
   ComponentCLKPtr comp = clkp->fcomp;
   EdgePtr ep1;
   while (comp)
   {
      if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
      {
         ep1 = comp->first_edge;
         for (int i = 0; i < comp->nedges; i++)
         {
            clkp->alt_lattice [lat(ep1->start)] = ep1->locpool + 1;
            ep1->frozen = false;
            ep1 = ep1->next;
         }
      }
      comp = comp->next;
   }

   // now go thru the edges looking for antiparallel edges at distance 1
   comp = clkp->fcomp;
   int kount = 0;
   // bfacf_lattice_info (clkp);

   while (comp)
   {
      if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
      {
         ep1 = comp->first_edge;
         for (int i = 0; i < comp->nedges; i++)
         {
            for (int dir = 0; dir < 4; dir++)
            {
               ivector test_loc;
               add_ivector(test_loc, ep1->start, increment_NEWSUD [turn [ep1->dir][dir]]);
               int IDp1 = clkp->alt_lattice [lat(test_loc)];
               if (IDp1)
               {
                  if (IDp1 < 1 || IDp1 >= clkp->poolsize)
                  {
                     printf("bad value: %d\n", IDp1);
                     exit(22);
                  }
                  EdgePtr ep2;
                  bool ok = true;
                  if (clk_recombo_direct)
                  {
                     ep2 = clkp->edgepool [IDp1 - 1]->prev;
                     if (clk_recombo_limit && ep1->comp == ep2->comp && clkp->ncomps > 1)
                        ok = false;
                     if (ok && anti [ep1->dir][ep2->dir] && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2) >= clk_min_arc_length_distance &&
                             !(ep1->next == ep2->prev || ep1->prev == ep2->next))
                     {
                        ++kount;
                        mark(ep1, ep2);
                     }
                  }
                  else
                  {
                     ep2 = clkp->edgepool [IDp1 - 1];
                     if (ep1->dir == ep2->dir && ep1->ID > ep2->ID && clk_arc_length_distance(ep1, ep2) >= clk_min_arc_length_distance &&
                             !(ep1->next == ep2->prev || ep1->prev == ep2->next))
                     {
                        ++kount;
                        mark(ep1, ep2);
                     }
                  }
               }
            }
            ep1 = ep1->next;
         }
      }
      comp = comp->next;
   }

   comp = clkp->fcomp;
   while (comp)
   {
      if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
      {
         ep1 = comp->first_edge;
         for (int i = 0; i < comp->nedges; i++)
         {
            clkp->alt_lattice [lat(ep1->start)] = 0;
            ep1 = ep1->next;
         }
      }
      comp = comp->next;
   }

   //  validate (clkp, "end of A");

   return kount;
}


#ifdef OLD_STUFF

void recombo_mark_frozen NOT USED(EdgePtr ep1, EdgePtr ep2)
{
   ep1->frozen = ep2->frozen = true;
}

int bfacf_perform_recombination NOT USED(CubicLatticeKnotPtr clkp)
{
   //recombo_index = 0;
   //recombo_clkp = clkp;
   //  validate (clkp, "start of bfacf_perform_recombination()");

   // find all the edges involved in recombo sites and mark them as "frozen"
   int numb = bfacf_perform_recombination(clkp, recombo_mark_frozen);
   if (numb < 1) return 0; // no recombo sites

   // the following block of code simply moves all the "frozen" edges to start of edgepool
   ComponentCLKPtr comp = clkp->fcomp;
   EdgePtr ep;
   clkp->nfrozen = 0;
   while (comp)
   {
      if (!(comp->flags & COMPONENT_CLK_FLAG_OPEN))
      {
         ep = comp->first_edge;
         for (int i = 0; i < comp->nedges; i++)
         {
            if (ep->frozen)
               freezeEdge(clkp, ep);
            ep->frozen = false; // so, the edge wasn't really frozen after all
            ep = ep->next;
         }
      }
      comp = comp->next;
   }

   int kount = 0;


   // we now have clkp->nfrozen edges at start of edge pool, choose one at random,
   // note that clkp->nfrozen may be less than 2 * numb because a given edge
   // might be in more than one recombo site
   ep = clkp->edgepool [rand_integer(0, clkp->nfrozen)];

   // move this edge to start of edge pool
   swapEdgesPool(clkp, ep->locpool, 0);

   //printf ("\n\n\n\n");

   // move all edges that are a distance 1 and antiparallel to ep 

   for (int i = 1; i < clkp->nfrozen; i++)
   {
      ivector incr;
      if (clk_recombo_direct)
      {
         sub_ivector(incr, ep->next->start, clkp->edgepool [i]->start);
         if (anti [ep->dir][clkp->edgepool [i]->dir] && clk_check_increment(incr))
         {
            swapEdgesPool(clkp, i, clkp->nfrozen + kount);
            //printf ("swap %d %d\n", i, clkp->nfrozen + kount);
            kount++;
         }
      }
      else
      {
         sub_ivector(incr, ep->start, clkp->edgepool [i]->start);
         if (ep->dir == clkp->edgepool [i]->dir && clk_check_increment(incr))
         {
            swapEdgesPool(clkp, i, clkp->nfrozen + kount);
            //printf ("swap %d %d\n", i, clkp->nfrozen + kount);
            kount++;
         }
      }
   }

   //showem (clkp, kount); clkp->nfrozen = 0; return 0;

   // kount edges were moved, and the candidate edges to recombo with ep
   // are numbered from clkp->nfrozen, ... clkp->nfrozen + kount - 1

   perform_recombination(clkp, ep, clkp->edgepool [rand_integer(clkp->nfrozen, clkp->nfrozen + kount)]);
   NOT USED

   clkp->nfrozen = 0;
   NOT USED
   return numb;
}

#endif

bool clk_verbose_knot_path_info(CubicLatticeKnotPtr clkp, char *filename)
{
   FILE *fp = fopen(filename, "w");
   if (!fp) return false;

   ComponentCLKPtr comp = clkp->fcomp;
   EdgePtr ep, epc;

   fprintf(fp, "ID, component, xstart, ystart, zstart, direction, locpool, prev, next, tclasp, parsite\n");

   while (comp)
   {
      ep = comp->first_edge;
      int tc, ps;
      for (int i = 0; i < comp->nedges; i++)
      {
         ivector loc;
         sub_ivector(loc, ep->start, clkp->loffset);
         tc = ps = -1;
         epc = is_tight_clasp(clkp, ep);
         if (epc) tc = epc->ID;
         epc = is_parsite(clkp, ep);
         if (epc) ps = epc->ID;

         fprintf(fp, "%d, %d, %d, %d, %d, %s, %d, %d, %d, %d, %d\n",
                 ep->ID, comp->ID,
                 loc [0], loc [1], loc [2],
                 clk_dir_name [ep->dir],
                 ep->locpool, ep->prev->ID, ep->next->ID,
                 tc, ps);
         ep = ep->next;
      }
      comp = comp->next;
   }

   fclose(fp);
   return true;
}

