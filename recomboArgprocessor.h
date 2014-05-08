/* 
 * File:   recomboArgprocessor.h
 * Author: kmo
 *
 * Created on February 28, 2013, 9:11 AM
 */

#ifndef RECOMBOARGPROCESSOR_H
#define	RECOMBOARGPROCESSOR_H

#include "argprocessor.h"
#include "clk.h"
#include "clkConformationAsList.h"

#include <string>

class recomboBaseArgprocessor : public ArgumentProcessor
{
protected:
   virtual bool processVector(std::vector<std::string> &arg);

   virtual void addOrderedArguments() = 0;
   
public:
   recomboBaseArgprocessor(const std::string& name);
   virtual ~recomboBaseArgprocessor();

   virtual const std::string& getOutputFileName() const;
   virtual const std::string& getLogFileName() const;
   virtual const std::string& getStatsFileName() const;

   virtual int getSeed() const;
   virtual bool seedWasSetUsingSystemTimer() const;

   virtual int getMinarc() const;
   virtual int getMaxarc() const;

private:
   // don't want copy constructor called
   recomboBaseArgprocessor(const recomboBaseArgprocessor& orig);
   // don't want default constructor called
   recomboBaseArgprocessor();

   std::string logFile;
   std::string statsFile;

   int s;
   bool setUsingSystemTimer;
   int min;
   int max;
};

class rerecomboArgprocessor: public recomboBaseArgprocessor
{
protected:
   virtual bool processVector(std::vector<std::string> &arg);

   virtual void addOrderedArguments();
   
public:
   rerecomboArgprocessor();
   virtual ~rerecomboArgprocessor();

   virtual const std::string& getAfterFileName() const;
   virtual const std::string& getBeforeFileName() const;

   virtual bool twocomonentlinks() const;

private:
   // don't want copy constructor called
   rerecomboArgprocessor(const rerecomboArgprocessor& orig);
   
   std::string afterFile;
   std::string beforeFile;
   bool processtwocomponentlink;

};

/**
 * Class to process arguments passed to main function of recombination 
 * application. 
 */
class recomboArgprocessor : public recomboBaseArgprocessor
{
protected:
   virtual bool processVector(std::vector<std::string> &arg);

   virtual void addOrderedArguments();

public:
   recomboArgprocessor();
   virtual ~recomboArgprocessor();

   virtual const std::string& getInitialConformationFileName() const;
   virtual const std::string& getAfterFileName() const;
   virtual const std::string& getBeforeFileName() const;

   //override
//      virtual const std::string& getOutputFileName() const;

   
   virtual double getZ() const;

   virtual long getNumberOfConformations() const;
   virtual long getNumberOfIterations() const;
   virtual bool useNumberOfIterations() const;
   virtual int getWarmup() const;
   virtual int getSampleRate() const;

   virtual bool hProtocol() const;
   virtual bool aProtocol() const;
   virtual bool bProtocol() const;

private:
   std::string afterFile;
   std::string beforeFile;

   double z;

   int w;
   int c;
   long nori;
   bool useN;

   bool h;
   bool a;
   bool b;

   // don't want copy constructor called
   recomboArgprocessor(const recomboArgprocessor& orig);

};

#endif	/* RECOMBOARGPROCESSOR_H */

