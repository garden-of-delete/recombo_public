/* 
 * File:   recomboArgprocessor.cpp
 * Author: kmo
 * 
 * Created on February 28, 2013, 9:11 AM
 */

#include "recomboArgprocessor.h"

#include <iostream>
#include <fstream>
#include <time.h>


using namespace std;

recomboBaseArgprocessor::recomboBaseArgprocessor(const string& name) : ArgumentProcessor(name)
{
   setUsingSystemTimer = true;
   time_t the_time;
   time(&the_time);
   s = the_time;

   // default settings
   min = 71;
   max = 77;

   //   addOrderedArguments();

   addOption("s", "SEED", "seed for the random number generator");
   addOption("arc", "ARC", "will insure that the arclength between recombination edges is exactly length ARC. If this option is set, then MINARC and MAXARC will be ignored. In general, only odd arclength of at least three make sense."); // Default is " + tostr(min) + ".");
   addOption("minarc", "MINARC", "will insure that the arclength between recombination edges is at least MINARC. If this option is set, then MAXARC must also be set. Default is " + tostr(max) + ".");
   addOption("maxarc", "MAXARC", "will insure that the arclength between recombination edges is at most MAXARC. If this option is set, then MINARC must also be set. Default is " + tostr(min) + ".");
}

recomboBaseArgprocessor::~recomboBaseArgprocessor() { }

bool recomboBaseArgprocessor::processVector(std::vector<std::string> &arg)
{
   if (!ArgumentProcessor::processVector(arg)) return false;

   if (!ordered[0].found) return false;

   logFile = getOutputFileName() + ".log";
   statsFile = getOutputFileName() + "_stats.txt";

   // was the seed specified, or should we use the system timer
   if (options["s"].found)
   {
      s = atoi(options["s"].content.c_str());
      setUsingSystemTimer = false;
   }

   // set the minimum and maximum arclengths
   if (options["arc"].found)
   {
      // both the min and max set to same
      min = max = atoi(options["arc"].content.c_str());
   }
   else if (options["minarc"].found || options["maxarc"].found)
   {
      // both minarc and maxarc must be set if either is set
      if (!options["minarc"].found || !options["maxarc"].found) return false;

      min = atoi(options["minarc"].content.c_str());
      max = atoi(options["maxarc"].content.c_str());
   }

   if (getMinarc() < 3)
   {
      cerr << "A minimum arclength of less than three does not make sense." << endl;
      return false;
   }
   if (getMaxarc() < getMinarc())
   {
      cerr << "A maximum arclength of less than the minimum arclength does not make sense." << endl;
      return false;
   }

   return true;
}

const string& recomboBaseArgprocessor::getOutputFileName() const
{
   if (ordered[1].found) return ordered[1].content;
   return ordered[0].content;
}

const string& recomboBaseArgprocessor::getLogFileName() const
{
   return logFile;
}

const string& recomboBaseArgprocessor::getStatsFileName() const
{
   return statsFile;
}

int recomboBaseArgprocessor::getSeed() const
{
   return s;
}

bool recomboBaseArgprocessor::seedWasSetUsingSystemTimer() const
{
   return setUsingSystemTimer;
}

int recomboBaseArgprocessor::getMinarc() const
{
   return min;
}

int recomboBaseArgprocessor::getMaxarc() const
{
   return max;
}

rerecomboArgprocessor::rerecomboArgprocessor() : recomboBaseArgprocessor("rerecombo")
{
   addOrderedArguments();

   addFlag("2", "this flag indicates that input file should be processed as two component links. If not set, then file will be processed as knots.");

}

rerecomboArgprocessor::~rerecomboArgprocessor() { }

bool rerecomboArgprocessor::processVector(vector<string> &arg)
{
   if (!recomboBaseArgprocessor::processVector(arg)) return false;

   if (!ordered[1].found) return false;
   beforeFile = ordered[0].content;
   afterFile = ordered[1].content + ".b";

   processtwocomponentlink = flags["2"].found;

   return true;
}

void rerecomboArgprocessor::addOrderedArguments()
{
   addOrdered("input", "filename of file containing sequence of conformations to recombine. This field is required.");
   addOrdered("output", "file name without extension for output files. After recombination conformations will be written to output.b and log will be written to output.log. This field is required.");
}

const string& rerecomboArgprocessor::getAfterFileName() const
{
   return afterFile;
}

const string& rerecomboArgprocessor::getBeforeFileName() const
{
   return beforeFile;
}

bool rerecomboArgprocessor::twocomonentlinks() const
{
   return processtwocomponentlink;
}

recomboArgprocessor::recomboArgprocessor() : recomboBaseArgprocessor("recombo")
{
   // default settings
   z = 0.208150;
   w = 1000000;
   nori = 60000000;
   useN = false;
   c = 25000;

   addOrderedArguments();

   addOption("z", "z-value", "z-value, default is " + tostr(z) + ".");
   addOption("c", "SAMPLE", "sampling rate, the minimum number of BFACF steps between samples. Default is " + tostr(c) + ".");
   addOption("w", "WARMUP", "do WARMUP iterations as a warmup. Default is " + tostr(w) + ".");
   addOption("i", "NITER", "do NITER iterations after warmup, and then stop. Default is " + tostr(nori) + ".");
   addOption("n", "NUMKNOTS", "ignore ITER and run simulation until NUMKNOTS are sampled.");

   addFlag("H", "use the simplest protocol. Filter sequence by recombination criterion without restarting. This is the default.");
   addFlag("A", "use the Metropolis-Hastings protocol. Take c steps, and test for recombination criterion. If new conformation does not satisfy criterion, restore old conformation.");
   addFlag("B", "use the windowed filter protocol.");

}

recomboArgprocessor::~recomboArgprocessor() { }

bool recomboArgprocessor::processVector(std::vector<std::string> &arg)
{
   if (!recomboBaseArgprocessor::processVector(arg)) return false;

   if (!ordered[0].found) return false;

   afterFile = getOutputFileName() + "_after.b";
   beforeFile = getOutputFileName() + "_before.b";

   // set z-value
   if (options["z"].found)
      z = atof(options["z"].content.c_str());

   // how many iterations for warmup?
   if (options["w"].found)
      w = atoi(options["w"].content.c_str());

   // do we want to do a fixed number of iterations, or run until we find a
   // number of conformations
   if (options["n"].found)
   {
      nori = atol(options["n"].content.c_str());
      useN = true;
   }
   else if (options["i"].found)
   {
      nori = atol(options["i"].content.c_str());
      useN = false;
   }

   // what sampling frequency?
   if (options["c"].found)
      c = atoi(options["c"].content.c_str());

   // check which protocol we want, from most conservative to most radical
   if (flags["A"].found)
   {
      h = false;
      a = true;
      b = false;
   }
   else if (flags["H"].found)
   {
      h = true;
      a = false;
      b = false;
   }
   else if (flags["B"].found)
   {
      h = false;
      a = false;
      b = true;
   }
   else
   {
      // default is to use the most conservative protocol
      h = false;
      a = true;
      b = false;
   }

   // do some sanity checking on the arguments
   if (useNumberOfIterations())
   {
      if (getNumberOfConformations() < 1)
      {
         cerr << "-n " << getNumberOfConformations() << " does not make sense. Must record at least one conformation." << endl;
         return false;
      }
   }
   else
   {
      if (getNumberOfIterations() < 1)
      {
         cerr << "-i " << getNumberOfIterations() << " does not make sense. Must do at least one iteration." << endl;
         return false;
      }
   }

   if (getSampleRate() < 1)
   {
      cerr << "-c " << getSampleRate() << " does not make sense. Must have a positive sampling rate." << endl;
      return false;
   }
   if (getWarmup() < 0)
   {
      cerr << "-w " << getWarmup() << " does not make sense. Must have a non-negative warmup." << endl;
      return false;
   }

   return true;
}

void recomboArgprocessor::addOrderedArguments()
{
   addOrdered("init", "filename of file containing initial conformation. This field is required.");
   addOrdered("output", "file name for output files, will create files output_before.b, output_after.b, output.log and output_stats.txt. If not specified, the same file name as init will be used.");
}

const string& recomboArgprocessor::getInitialConformationFileName() const
{
   return ordered[0].content;
}

const string& recomboArgprocessor::getAfterFileName() const
{
   return afterFile;
}

const string& recomboArgprocessor::getBeforeFileName() const
{
   return beforeFile;
}

double recomboArgprocessor::getZ() const
{
   return z;
}

long recomboArgprocessor::getNumberOfConformations() const
{
   return nori;
}

long recomboArgprocessor::getNumberOfIterations() const
{
   return nori;
}

bool recomboArgprocessor::useNumberOfIterations() const
{
   return useN;
}

int recomboArgprocessor::getWarmup() const
{
   return w;
}

int recomboArgprocessor::getSampleRate() const
{
   return c;
}

bool recomboArgprocessor::hProtocol() const
{
   return h;
}

bool recomboArgprocessor::aProtocol() const
{
   return a;
}

bool recomboArgprocessor::bProtocol() const
{
   return b;
}

