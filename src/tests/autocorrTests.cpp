#include "testFramework.h"

#include "autocorr.h"

#include <cmath>

bool autocorrelation(double &mean, double &error_mean,
      double &autocorr, double &error_autocorr,
      int &window_size,
      const double *data, int num_values, bool verbose)
{

    double e, w, z;
    double s, sx, sy, sxy, sxx;
    int k, l, i;
    double a, t, sa, sb;
    double mean_sq;

    // Calculate the mean of the stochastic variable
    mean = 0.0;
    for (i = 0; i < num_values; i++)
        mean += data [i];
    mean = mean / (double) num_values;

    // Calculate the mean of the square stochastic variable
    mean_sq = 0.0;
    for (i = 0; i < num_values; i++)
        mean_sq += data [i] * data [i];
    mean_sq = mean_sq / (double) num_values;

    // Calculate the autocorrelation function now
    mean_sq -= mean * mean;
    e = 0.;
    s = 0.;
    sx = 0.;
    sy = 0.;
    sxy = 0.;
    sxx = 0.;

    for (k = 0; k < num_values; k++)
    {
        w = 0.;
        z = 0.;
        for (l = 0; l < num_values - 1 - k; l++)
        {
            w += 1.0;
            z += 2.0 * (data [l] - mean) * (data [l + k + 1] - mean);
        }
        z = z / w;
        if (z > 0.)
            e += z;
        else
            break;

        a = z / (2.0 * mean_sq); /* This is the value of the autocorrelation function.  If we in its
				 then we can try to estimate its exponential tail. */

        // RGS: the logic of the following is a bit sketchy

        s += 1.0;
        sx += (double) k + 1.0;
        sxx += ((double) k + 1.0) * ((double) k + 1.0);
        sy += log(z / 2.0);
        sxy += ((double) k + 1.0) * log(z / 2.0);
        if (s > 3.)
        { // RGS: why 3?
            sa = (s * sxy - sx * sy) / (s * sxx - sx * sx);
            sa = -1.0 / sa;
            sb = (sxx * sy - sx * sxy) / (s * sxx - sx * sx);
        } else
        {
            sa = 0.;
            sb = 0.;
        }

        if (s == 30)
        { // RGS: why 30?
            s = 0.;
            sx = 0.;
            sxx = 0.;
            sy = 0.;
            sxy = 0.;
        }

        t = (mean_sq + e) / (2.0 * mean_sq);

//        if (verbose)
//            cout << "I=" << k << ", C=" << a << ", T=" << t << ", E=" << sa << ", K="
//              << sb << endl;

        if (k > (int) (1000.0 * t))
        {
            if (verbose)
                cerr << "The window is " << 1000 << " times larger than the autocorrelation time." << endl;
            return false;
        }
    }

    window_size = k;
    autocorr = t = (mean_sq + e) / (2.0 * mean_sq);
    error_autocorr = sb = t * sqrt(2.0 * (2.0 * (double) k + 1.0) / (double) num_values);
    error_mean = sqrt(2.0 * mean_sq * t / (double) num_values);
    return true;
}


const double autocorrTestDataNumVals = 30;

const double autocorrTestData[] = {
    1.2190079223,
    -8.9262668602,
    7.7765452396,
    -1.3448207639,
    -7.7495624404,
    -1.7062141187,
    -7.2927681822,
    3.8728410844,
    -6.2263957225,
    -6.1935190763,
    7.8775635548,
    4.2208929267,
    -1.3419686165,
    8.3451695926,
    6.3461351953,
    3.8499901164,
    7.0610383339,
    -9.9842282571,
    0.0627100654,
    -8.4068057034,
    -1.7178643402,
    -4.5115231071,
    -4.5541757066,
    0.8381974231,
    2.4940358102,
    -9.5345366467,
    -1.4188661613,
    -1.5643732436,
    2.2443082184,
    0.4945876356
};

#define EPSILON 1.0e-10
#define DELTA 1.0e-4

bool testAutocorr()
{
    // modified form newsimpletest autocorr test.
    vector<double> data;
    for (int i = 0; i < autocorrTestDataNumVals; i++)
        data.push_back(autocorrTestData[i]);
    // Expect 30 data values.
    ASSERT_MESSAGE("number of data values incorrect", data.size() == autocorrTestDataNumVals);
    autocorr ac;
    double mean = ac.computeMean(autocorrTestData, autocorrTestDataNumVals);
    double meanByVector = ac.computeMean(data);
    double expectedMean = -0.8590288609;
    
    ASSERT_DOUBLES_EQUAL_MESSAGE("mean incorrect", mean, expectedMean, EPSILON);
    ASSERT_DOUBLES_EQUAL_MESSAGE("means should be same computed by both methods", mean, meanByVector, EPSILON);

    double mean_sq = ac.computeMeanSquare(autocorrTestData, autocorrTestDataNumVals);
    double mean_sqByVector = ac.computeMeanSquare(data);
    double expectedMean_sq = 30.9991617759;
    double expectedVariance = 30.261231192;
    double var = ac.computeVariance(autocorrTestData, autocorrTestDataNumVals);
    double varByVector = ac.computeVariance(data);

    ASSERT_DOUBLES_EQUAL_MESSAGE("mean square incorrect",mean_sq, expectedMean_sq, EPSILON);
    ASSERT_DOUBLES_EQUAL_MESSAGE("mean squares should be same computed by both methods", mean_sq, mean_sqByVector, EPSILON);
    ASSERT_DOUBLES_EQUAL_MESSAGE("variance incorrect", var, expectedVariance, EPSILON);
    ASSERT_DOUBLES_EQUAL_MESSAGE("variance should be same computed by both methods", var, varByVector, EPSILON);

    // compute the values exactly as in Stu Whittington's program 
    double meanSW, mean_sqSW;
    int i, num_values;

    num_values = autocorrTestDataNumVals;

    // Calculate the mean of the stochastic variable
    meanSW = 0.0;
    for (i = 0; i < num_values; i++)
        meanSW += data [i];
    meanSW = meanSW / (double) num_values;

    // Calculate the mean of the square stochastic variable
    mean_sqSW = 0.0;
    for (i = 0; i < num_values; i++)
        mean_sqSW += data [i] * data [i];
    mean_sqSW = mean_sqSW / (double) num_values;

    ASSERT_DOUBLES_EQUAL_MESSAGE("means should be same computed by SW", mean, meanSW, EPSILON);
    ASSERT_DOUBLES_EQUAL_MESSAGE("mean squares should be same computed by SW", mean_sq, mean_sqSW, EPSILON);

    // Calculate the autocorrelation function now
    mean_sqSW -= meanSW * meanSW;

    // mean_sqSW should be the variance now
    ASSERT_DOUBLES_EQUAL_MESSAGE("variance should be same computed by SW", var, mean_sqSW, EPSILON);

    autocorrInfo aci(1, 2, 3, 4, 5);
//    cout << aci << endl;
    autocorrInfo aci2;
    aci2 = aci;
//    cout << aci2 << endl;

    // check that default == operator works for autocorrInfo
    ASSERT_MESSAGE("Redefine ==, please.", aci == aci2);

    aci2.autocorr = 0.5;
    ASSERT_MESSAGE("Redefine ==, please.", aci != aci2);

    // test algorithm in class versus the algorithm in class
    aci.computationError = autocorrelation(aci.mean, aci.error_mean,
          aci.autocorr, aci.error_autocorr,
          aci.window_size,
          autocorrTestData, autocorrTestDataNumVals, false);

//    cout << aci << endl;

    aci2 = ac.autocorrelation(autocorrTestData, autocorrTestDataNumVals, false);

//    cout << aci2 << endl;

    ASSERT_MESSAGE("autocorr should produce same results as SW", aci == aci2);

    autocorrInfo aci3 = ac.autocorrelation(data, false);
    ASSERT_MESSAGE("autocorr should produce same results by either method", aci == aci3);
    
    return true;
}
