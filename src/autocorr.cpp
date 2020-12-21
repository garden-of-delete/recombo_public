/* 
 * File:   autocorr.cpp
 * Author: kmo
 * 
 * Created on September 21, 2012, 8:00 AM
 */

#include "autocorr.h"

#include <cmath>

autocorr::autocorr() : why3(3.0), why30(30), maxwindow(1000) { }

autocorr::autocorr(const autocorr& orig) : why3(orig.why3), why30(orig.why30), maxwindow(orig.maxwindow) { }

autocorr::~autocorr() { }

double autocorr::computeMean(const vector<double>& data) const
{
   double mean = 0.0;
   unsigned num_values = data.size();
   for (int i = 0; i < num_values; i++)
      mean += data [i];
   mean = mean / (double) num_values;
   return mean;
}

double autocorr::computeMean(const double* data, unsigned num_values) const
{
   double mean = 0.0;
   for (int i = 0; i < num_values; i++)
      mean += data [i];
   mean = mean / (double) num_values;
   return mean;
}

double autocorr::computeMeanSquare(const vector<double>& data) const
{
   double mean_sq = 0.0;
   unsigned num_values = data.size();
   for (int i = 0; i < num_values; i++)
      mean_sq += data [i] * data [i];
   mean_sq = mean_sq / (double) num_values;
   return mean_sq;
}

double autocorr::computeMeanSquare(const double* data, unsigned num_values) const
{
   double mean_sq = 0.0;
   for (int i = 0; i < num_values; i++)
      mean_sq += data [i] * data [i];
   mean_sq = mean_sq / (double) num_values;
   return mean_sq;
}

double autocorr::computeVariance(const vector<double>& data) const
{
   double mean = computeMean(data);
   return computeMeanSquare(data) - mean*mean;
}

double autocorr::computeVariance(const double* data, unsigned num_values) const
{
   double mean = computeMean(data, num_values);
   return computeMeanSquare(data, num_values) - mean*mean;
}

autocorrInfo autocorr::autocorrelation(const double *data, int num_values, bool verbose) const
{

   double e, w, z;
   double s, sx, sy, sxy, sxx;
   int k, l, i;
   double a, t, sa, sb;
   double mean = computeMean(data, num_values);
   autocorrInfo info;

   // Calculate the mean of the square stochastic variable
   double variance = computeMeanSquare(data, num_values) - mean * mean;
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

      a = z / (2.0 * variance); /* This is the value of the autocorrelation function.  If we in its
				 then we can try to estimate its exponential tail. */

      // RGS: the logic of the following is a bit sketchy

      s += 1.0;
      sx += (double) k + 1.0;
      sxx += ((double) k + 1.0) * ((double) k + 1.0);
      sy += log(z / 2.0);
      sxy += ((double) k + 1.0) * log(z / 2.0);
      if (s > why3)
      { // RGS: why 3?
         sa = (s * sxy - sx * sy) / (s * sxx - sx * sx);
         sa = -1.0 / sa;
         sb = (sxx * sy - sx * sxy) / (s * sxx - sx * sx);
      }
      else
      {
         sa = 0.;
         sb = 0.;
      }

      if (s == why30)
      { // RGS: why 30?
         s = 0.;
         sx = 0.;
         sxx = 0.;
         sy = 0.;
         sxy = 0.;
      }

      t = (variance + e) / (2.0 * variance);

      if (verbose)
         cout << "I=" << k << ", C=" << a << ", T=" << t << ", E=" << sa << ", K="
              << sb << endl;

      if (k > (int) (maxwindow * t))
      {
         if (verbose)
            cerr << "The window is " << maxwindow << " times larger than the autocorrelation time." << endl;
         info.computationError = true;
         return info;
      }
   }

   info.mean = mean;
   info.window_size = k;
   info.autocorr = t = (variance + e) / (2.0 * variance);
   info.error_autocorr = sb = t * sqrt(2.0 * (2.0 * (double) k + 1.0) / (double) num_values);
   info.error_mean = sqrt(2.0 * variance * t / (double) num_values);
//   cout << info << endl;
   return info;
}

autocorrInfo autocorr::autocorrelation(const vector<double>& data, bool verbose) const
{
   int num_values = data.size();
   double e, w, z;
   double s, sx, sy, sxy, sxx;
   int k, l, i;
   double a, t, sa, sb;
   double mean = computeMean(data);
   autocorrInfo info;

   // Calculate the mean of the square stochastic variable
   double variance = computeMeanSquare(data) - mean * mean;
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

      a = z / (2.0 * variance); /* This is the value of the autocorrelation function.  If we in its
				 then we can try to estimate its exponential tail. */

      // RGS: the logic of the following is a bit sketchy

      s += 1.0;
      sx += (double) k + 1.0;
      sxx += ((double) k + 1.0) * ((double) k + 1.0);
      sy += log(z / 2.0);
      sxy += ((double) k + 1.0) * log(z / 2.0);
      if (s > why3)
      { // RGS: why 3?
         sa = (s * sxy - sx * sy) / (s * sxx - sx * sx);
         sa = -1.0 / sa;
         sb = (sxx * sy - sx * sxy) / (s * sxx - sx * sx);
      }
      else
      {
         sa = 0.;
         sb = 0.;
      }

      if (s == why30)
      { // RGS: why 30?
         s = 0.;
         sx = 0.;
         sxx = 0.;
         sy = 0.;
         sxy = 0.;
      }

      t = (variance + e) / (2.0 * variance);

      if (verbose)
         cout << "I=" << k << ", C=" << a << ", T=" << t << ", E=" << sa << ", K="
              << sb << endl;

      if (k > (int) (maxwindow * t))
      {
         if (verbose)
            cerr << "The window is " << maxwindow << " times larger than the autocorrelation time." << endl;
         info.computationError = true;
         return info;
      }
   }

   info.mean = mean;
   info.window_size = k;
   info.autocorr = t = (variance + e) / (2.0 * variance);
   info.error_autocorr = sb = t * sqrt(2.0 * (2.0 * (double) k + 1.0) / (double) num_values);
   info.error_mean = sqrt(2.0 * variance * t / (double) num_values);
//   cout << info << endl;
   return info;
}

autocorrInfo::autocorrInfo() : computationError(false) { }

autocorrInfo::autocorrInfo(double autocorr, double error_autocorr, double mean, double error_mean, int window_size) :
mean(mean), error_mean(error_mean), autocorr(autocorr),
error_autocorr(error_autocorr), window_size(window_size),
computationError(false) { }

autocorrInfo::autocorrInfo(const autocorrInfo & orig) :
mean(orig.mean), error_mean(orig.error_mean), autocorr(orig.autocorr),
error_autocorr(orig.error_autocorr), window_size(orig.window_size),
computationError(orig.computationError) { }

autocorrInfo::~autocorrInfo() { }

bool autocorrInfo::operator==(const autocorrInfo& aci) const
{
   return autocorr == aci.autocorr && error_autocorr == aci.error_autocorr
           && mean == aci.mean && error_mean == aci.error_mean
           && window_size == aci.window_size;
}

bool autocorrInfo::operator!=(const autocorrInfo& aci) const
{
   return !(*this == aci);
}

autocorrInfo& autocorrInfo::operator=(const autocorrInfo& aci)
{
   autocorr = aci.autocorr;
   error_autocorr = aci.error_autocorr;
   mean = aci.mean;
   error_mean = aci.error_mean;
   window_size = aci.window_size;
   computationError = aci.computationError;
   return *this;
}

ostream &operator<<(ostream& os, const autocorrInfo& aci)
{
   return os << aci.autocorr << " " << aci.error_autocorr << " "
           << aci.mean << " " << aci.error_mean;
}

