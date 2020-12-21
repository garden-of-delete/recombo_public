/* 
 * File:   autocorr.h
 * Author: kmo
 *
 * Created on September 21, 2012, 8:00 AM
 */

#ifndef AUTOCORR_H
#define	AUTOCORR_H

#include <iostream>
#include <vector>

using namespace std;

struct autocorrInfo;

/**
 * A C++ class using modified by Reuben Brasher from code modified by Rob 
 * Scharein from the original C program by Buks van Rensburg and Enzo Orlandini
 * Thanks to Stu Wittington for giving access to the original version.
 *
 * </p>
 * Computes windowed autocorrelation and mean of a stochastic variable.
 */
class autocorr
{
public:
   /**
    * Default constructor for autocorr class. Sets default value for max window 
    * size to be 1000.
    */
   autocorr();

   /**
    * Copy constructor for autocorr class. Copies max window size.
    * @param orig Object to copy.
    */
   autocorr(const autocorr& orig);

   /**
    * Destructor for autocorr class.
    */
   virtual ~autocorr();

   /**
    * Computes the mean of a stochastic variable.
    * @param data Values of a stochastic variable.
    * @return Mean of data.
    */
   double computeMean(const vector<double>& data) const;

   /**
    * Computes the mean of a stochastic variable.
    * @param data Values of a stochastic variable.
    * @param num_values Length of data array.
    * @return Mean of data.
    */
   double computeMean(const double* data, unsigned num_values) const;

   /**
    * Computes the mean of square of stochastic variable.
    * @param data Values of a stochastic variable.
    * @return Mean of squares of data.
    */
   double computeMeanSquare(const vector<double>& data) const;

   /**
    * Computes the mean of square of stochastic variable.
    * @param data Values of a stochastic variable.
    * @param num_values Length of data array.
    * @return Mean of squares of data.
    */
   double computeMeanSquare(const double* data, unsigned num_values) const;

   /**
    * Computes the variance of a stochastic variable normalized by one over the
    * number of values. This is the standard biased variance.
    * @param data Values of a stochastic variable.
    * @return Variance of data.
    */
   double computeVariance(const vector<double>& data) const;

   /**
    * Computes the variance of a stochastic variable normalized by one over the
    * number of values. This is the standard biased variance.
    * @param data Values of a stochastic variable.
    * @param num_values Length of data array.
    * @return Variance of data.
    */
   double computeVariance(const double* data, unsigned num_values) const;

   /**
    * Compute autocorrelation time and average of a stochastic variable.
    * @param data Values of a stochastic variable.
    * @param num_values Length of data array.
    * @param verbose true for long messages, false otherwise.
    * @return Computed integrated autocorrelation time, standard deviation of 
    * autocorrelation time, mean of stochastic variable, and standard 
    * deviation of mean.
    */
   autocorrInfo autocorrelation(const double *data, int num_values, bool verbose) const;

   /**
    * Compute autocorrelation time and average of a stochastic variable.
    * @param data Values of a stochastic variable.
    * @param verbose true for long messages, false otherwise.
    * @return Computed integrated autocorrelation time, standard deviation of 
    * autocorrelation time, mean of stochastic variable, and standard 
    * deviation of mean.
    */
   autocorrInfo autocorrelation(const vector<double>& data, bool verbose) const;

private:
   // these are mystery variables for which RGS asked why
   double why3;
   double why30;
   int maxwindow;
};

struct autocorrInfo
{
   double autocorr;
   double error_autocorr;
   double mean;
   double error_mean;
   int window_size;
   // TODO: remove computationError field and replace with throwing an exception
   bool computationError;

   autocorrInfo();
   autocorrInfo(double autocorr, double error_autocorr, double mean, double error_mean, int window_size);
   autocorrInfo(const autocorrInfo & orig);
   virtual ~autocorrInfo();

   bool operator==(const autocorrInfo & aci) const;
   bool operator!=(const autocorrInfo & aci) const;
   autocorrInfo& operator=(const autocorrInfo & aci);

   friend ostream &operator<<(ostream&, const autocorrInfo&);
};

#endif	/* AUTOCORR_H */

