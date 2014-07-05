/* 
 * File:   runRecombo.h
 * Author: kmo
 *
 * Created on February 28, 2013, 9:11 AM
 */

#ifndef RUNRECOMBO_H
#define	RUNRECOMBO_H

#include <iostream>
#include <fstream>

#include "recomboArgprocessor.h"
#include "clkConformationBfacf3.h"
#include "pseudorandom.h"

/**
 * A class to handle the details of recombination experiments.
 */
class runRecombo
{
protected:
   recomboArgprocessor& args;

   std::fstream logfile;
   std::fstream statsfile;
   std::fstream before;
   std::fstream after;

   clkConformationBfacf3 *conformation;

   bool finished(long conformationsRecorded, long iterationsDone) const;
   void warmup();

   time_t starttime;

   void timestamp(std::ostream& os, time_t time) const;

   pseudorandom siteSelector;

   void recordConformation(std::ostream& os, clkConformationBfacf3 *conformation);

   clkConformationAsList *savedComponent0, *savedComponent1;
   void markConformation();
   void restoreConformation();

   //   clkConformationAsList initialCoords;
   int numberOfComponents;
   clkConformationAsList initialComp0, initialComp1;

public:
   runRecombo(recomboArgprocessor& args);
   virtual ~runRecombo();

   void writeLog(std::ostream& os) const;
   bool openfiles();
   void closefiles();

   virtual void dorecombo() = 0;

   bool readInitialCoords();
   bool readInitialCoords(std::istream& is);

   int getNumberOfComponents() const;
   const clk& getInitialComp0() const;
   const clk& getInitialComp1() const;

private:
   // don't want copy constructor called
   runRecombo(const runRecombo& orig);
};

class runRecomboA : public runRecombo
{
public:
   runRecomboA(recomboArgprocessor& args);
   virtual ~runRecomboA();

   virtual void dorecombo();
};

class runRecomboH : public runRecombo
{
public:
   runRecomboH(recomboArgprocessor& args);
   virtual ~runRecomboH();

   virtual void dorecombo();
};

class runRecomboB : public runRecombo
{
public:
   runRecomboB(recomboArgprocessor& args);
   virtual ~runRecomboB();

   virtual void dorecombo();
};

class runRerecombo
{
protected:
   rerecomboArgprocessor& args;

   std::fstream logfile;
   std::fstream statsfile;
   std::fstream before;
   std::fstream after;

   time_t starttime;
   void timestamp(std::ostream& os, time_t time) const;

   pseudorandom siteSelector;

   void recordConformation(std::ostream& os, clkConformationBfacf3 *conformation);
   clkConformationBfacf3* readConformation(std::istream& is);

public:
   runRerecombo(rerecomboArgprocessor& args);
   virtual ~runRerecombo();

   void writeLog(std::ostream& os) const;
   bool openfiles();
   void closefiles();

   virtual void dorecombo();

};

#endif	/* RUNRECOMBO_H */

